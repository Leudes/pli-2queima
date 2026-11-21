"""
Modelo de Programação Linear Inteira (PLI) para o Número de 2-Queima
(2-Burning Number, b2(G)) via variável z que codifica o instante final
usando a variável indicadora f_t, com WARM START vindo de um algoritmo guloso.

Variáveis:
  x[v,t]  -> 1 se o vértice v está queimado até o instante t (t = 0,...,T)
  s[v,t]  -> 1 se v é a fonte de queima escolhida no instante t (t = 1,...,T)
  f[t]    -> 1 se o processo termina exatamente no instante t
  z       -> instante final (valor objetivo)

Requer: JuMP, HiGHS, Graphs (todos instaláveis via Pkg).
"""

using JuMP
using Graphs
using HiGHS

# ---------------------------------------------------------------------
# Valores padrão. Todos podem ser alterados pela linha de comando.
# ---------------------------------------------------------------------
const TIME_LIMIT     = 300.0  # limite de tempo do solver por grafo (segundos)
const WARM_START     = true   # fornece a solução gulosa como solução inicial ao solver
const SHRINK_HORIZON = true   # usa T = valor guloso (limitante superior válido) em vez de T = n
const SILENT         = true   # false mostra o log do solver (útil p/ conferir o MIP start)
const SOLVER_THREADS = 1      # evita que execuções disputem todos os núcleos
const MAX_BINARY_VARIABLES = 1_000_000
const CABECALHO_CSV =
    "graph,vertices,edges,greedy_method,greedy_burning_number," *
    "pli_burning_number,objective_bound,gap,time_greedy,time_model," *
    "time_solver,time_total,#sources,status"


# ---------------------------------------------------------------------
# Função que lê um grafo simples de um arquivo texto no formato de lista de arestas 
# ---------------------------------------------------------------------
function read_simple_graph(filepath::String)::SimpleGraph
    unique_edges = Set{Tuple{Int,Int}}()
    vertex_labels = Set{Int}()
    declared_vertices = nothing

    for line in eachline(filepath)
        line = strip(line)
        (isempty(line) || startswith(line, '#') || startswith(line, '%')) && continue
        parts = split(line)

        # Também aceita o cabeçalho DIMACS "p edge n m".
        if parts[1] == "p"
            length(parts) >= 4 || error("Cabeçalho DIMACS inválido em $filepath")
            declared_vertices = parse(Int, parts[end - 1])
            continue
        end

        # Também aceita linhas DIMACS "e u v".
        offset = parts[1] == "e" ? 1 : 0
        length(parts) >= offset + 2 || continue
        u = parse(Int, parts[offset + 1])
        v = parse(Int, parts[offset + 2])
        push!(vertex_labels, u)
        push!(vertex_labels, v)
        u != v && push!(unique_edges, u < v ? (u, v) : (v, u))
    end

    if declared_vertices !== nothing
        # Arquivos DIMACS são indexados de 1 a n.
        union!(vertex_labels, 1:declared_vertices)
    end
    isempty(vertex_labels) && error("O arquivo não possui vértices: $filepath")

    all_vertices = sort!(collect(vertex_labels))
    vertex_map = Dict(v => i for (i, v) in enumerate(all_vertices))

    g = SimpleGraph(length(all_vertices))
    for (u, v) in unique_edges
        add_edge!(g, vertex_map[u], vertex_map[v])
    end
    return g
end


# ---------------------------------------------------------------------
# Algoritmo guloso para o 2-queima
# ---------------------------------------------------------------------
#
# Regra do processo (a mesma codificada nas restrições R1-R8):
#   * Rodada t: primeiro ocorre a propagação a partir do estado em t-1
#     (um vértice não queimado queima se tiver >= 2 vizinhos queimados);
#     depois, se ainda restar vértice não queimado, escolhe-se uma nova
#     fonte entre os não queimados.
#   * O processo termina na primeira rodada t em que todos estão queimados.
#
# A escolha da fonte é a mesma da heurística inicial oficial do GVNS:
# combina grau relativo e quantidade de gatilhos criados.

"Um passo de propagação: retorna o novo estado de queima a partir de `burned`."
function propagate_step(g::AbstractGraph, burned::AbstractVector{Bool})
    new = copy(burned)
    for v in vertices(g)
        burned[v] && continue
        c = 0
        for u in neighbors(g, v)
            if burned[u]
                c += 1
                c >= 2 && break
            end
        end
        c >= 2 && (new[v] = true)
    end
    return new
end

"""
    greedy_2burning(g)

Guloso usado como solução inicial em todos os grafos. Em cada rodada, escolhe o vértice
que maximiza o maior valor entre:

  * grau relativo: quantidade de vizinhos que ainda não queimou;
  * gatilhos: vizinhos que passariam de um para dois vizinhos queimados.

Cada rodada percorre os vértices e arestas do grafo. A regra é a mesma
empregada na solução inicial oficial do GVNS.
"""
function greedy_2burning(g::AbstractGraph)
    n = nv(g)
    burned = falses(n)
    states = [copy(burned)]
    sources = Tuple{Int,Int}[]
    t = 0

    while true
        t += 1
        new = propagate_step(g, burned)

        if all(new)
            push!(states, copy(new))
            break
        end

        burned_neighbors = zeros(Int, n)
        for v in vertices(g)
            new[v] && continue
            burned_neighbors[v] = count(u -> new[u], neighbors(g, v))
        end

        best = 0
        best_score = (-1, -1, -1)
        for v in vertices(g)
            new[v] && continue

            relative_degree = degree(g, v) - burned_neighbors[v]
            triggers = count(
                u -> !new[u] && burned_neighbors[u] == 1,
                neighbors(g, v),
            )
            score = (max(relative_degree, triggers), triggers, relative_degree)

            if score > best_score
                best = v
                best_score = score
            end
        end

        best != 0 || error("Não foi possível escolher uma nova fonte")
        new[best] = true
        push!(sources, (t, best))
        push!(states, copy(new))
        burned = new
        all(new) && break
    end

    return t, sources, states
end


# ---------------------------------------------------------------------
# Construção do modelo
# ---------------------------------------------------------------------

"""
    build_2burning_model(g::AbstractGraph, T::Int;
                         silent::Bool=false,
                         time_limit::Union{Nothing,Real}=nothing)

Constrói o modelo PLI (R1-R11) para o grafo `g` com horizonte de tempo `T`.
`T` deve ser um limitante superior válido para o número de queima
(por exemplo, T = nv(g), ou o valor obtido pelo algoritmo guloso).
Se `time_limit` for informado, o solver para após esse número de segundos.

Retorna o `Model` do JuMP pronto para ser otimizado.
"""
function build_2burning_model(g::AbstractGraph, T::Int;
                              silent::Bool=false,
                              time_limit::Union{Nothing,Real}=nothing,
                              threads::Int=SOLVER_THREADS)
    V = collect(vertices(g))

    model = Model(HiGHS.Optimizer)
    silent && set_silent(model)
    time_limit !== nothing && set_time_limit_sec(model, Float64(time_limit))
    set_optimizer_attribute(model, "threads", threads)

    d(v) = degree(g, v)
    N(v) = neighbors(g, v)

    # ---------------- Variáveis ----------------
    @variable(model, x[v in V, t in 0:T], Bin)
    @variable(model, s[v in V, t in 1:T], Bin)
    @variable(model, f[t in 1:T], Bin)
    @variable(model, z >= 1, Int)

    # ---------------- R1: x_{v,0} = 0 ----------------
    @constraint(model, R1[v in V], x[v, 0] == 0)

    # ---------------- R2: x_{v,t} >= x_{v,t-1} (monotonicidade) ----------------
    @constraint(model, R2[v in V, t in 1:T], x[v, t] >= x[v, t-1])

    # ---------------- R3: x_{v,t} >= s_{v,t} ----------------
    @constraint(model, R3[v in V, t in 1:T], x[v, t] >= s[v, t])

    # ---------------- R4: cada vértice é fonte no máximo uma vez ----------------
    @constraint(model, R4[v in V], sum(s[v, t] for t in 1:T) <= 1)

    # ---------------- R5: no máximo uma fonte por instante ----------------
    @constraint(model, R5[t in 1:T], sum(s[v, t] for v in V) <= 1)

    # ---------------- R6: só pode ser fonte se ainda não queimado ----------------
    @constraint(model, R6[v in V, t in 1:T], s[v, t] <= 1 - x[v, t-1])

    # ---------------- R7: propagação (upper bound de queima em 2 saltos) --------
    @constraint(model, R7[v in V, t in 1:T],
        x[v, t] <= x[v, t-1] + s[v, t] + 0.5 * sum(x[u, t-1] for u in N(v)))

    # ---------------- R8: só p/ vértices com grau >= 2 ---------------------------
    @constraint(model, R8[v in V, t in 1:T; d(v) >= 2],
        0.5 + d(v) * x[v, t] >=
            0.5 * sum(x[u, t-1] for u in N(v)) - d(v) * (x[v, t-1] + s[v, t]))

    # ---------------- R9: exatamente um instante final ----------------
    @constraint(model, R9, sum(f[t] for t in 1:T) == 1)

    # ---------------- R10: z = instante final ----------------
    @constraint(model, R10, z == sum(t * f[t] for t in 1:T))

    # ---------------- R11: todo vértice deve estar queimado no instante final ---
    @constraint(model, R11[v in V, t in 1:T], x[v, t] >= f[t])

    # ---------------- Objetivo ----------------
    @objective(model, Min, z)

    return model
end

# ---------------------------------------------------------------------
# Warm start
# ---------------------------------------------------------------------

"""
    set_warm_start!(model, g, T, z_g, sources, states)

Fornece ao modelo a solução gulosa como solução inicial (MIP start),
preenchendo TODAS as variáveis (x, s, f, z). Requer T >= z_g.
Para t > z_g, o estado fica todo queimado (x = 1) e s = 0.
"""
function set_warm_start!(model, g::AbstractGraph, T::Int, z_g::Int,
                         sources::Vector{Tuple{Int,Int}},
                         states::Vector{<:AbstractVector{Bool}})
    T >= z_g || error("Horizonte T=$T menor que o valor guloso z=$z_g")

    x, s, f, z = model[:x], model[:s], model[:f], model[:z]

    for v in vertices(g), t in 0:T
        tt = min(t, z_g)
        set_start_value(x[v, t], states[tt+1][v] ? 1.0 : 0.0)
    end

    for v in vertices(g), t in 1:T
        set_start_value(s[v, t], 0.0)
    end
    for (t, v) in sources
        set_start_value(s[v, t], 1.0)
    end

    for t in 1:T
        set_start_value(f[t], t == z_g ? 1.0 : 0.0)
    end

    set_start_value(z, float(z_g))
    return nothing
end

# ---------------------------------------------------------------------
# Extração da solução
# ---------------------------------------------------------------------

"""
    extract_solution(model, g)

Após `optimize!(model)`, extrai:
  - b2 : valor de z (ótimo, ou melhor solução encontrada se houve limite de tempo)
  - sources : vetor de (t, v) com as fontes escolhidas em ordem de t
"""
function extract_solution(model, g::AbstractGraph)
    b2 = round(Int, value(model[:z]))

    s = model[:s]
    sources = Tuple{Int,Int}[]
    for t in 1:last(axes(s, 2))
        for v in vertices(g)
            if value(s[v, t]) > 0.5
                push!(sources, (t, v))
            end
        end
    end
    sort!(sources, by = first)

    return b2, sources
end

# ---------------------------------------------------------------------
# Saída em CSV
# ---------------------------------------------------------------------

"""
    init_csv(path)

Cria o arquivo CSV com o cabeçalho, caso ele ainda não exista.
Se já existir, as novas linhas são anexadas ao final.
"""
function init_csv(path::AbstractString)
    if !isfile(path)
        open(path, "w") do io
            println(io, CABECALHO_CSV)
        end
    end
end

"""
    append_csv_row(path, graph_name, burning_number, time, status, greedy_bn, sources)

Anexa uma linha ao CSV. `burning_number` pode ser `missing` quando o solver
não encontrou nenhuma solução viável (campo fica vazio). `greedy_bn` é o valor
obtido pelo algoritmo guloso.
"""
function append_csv_row(path::AbstractString, graph_name::AbstractString;
                        vertices::Int, edges::Int, greedy_method::AbstractString,
                        greedy_bn, burning_number, objective_bound_value,
                        gap_value, time_greedy::Real, time_model::Real,
                        time_solver::Real, time_total::Real, sources, status)
    bn = ismissing(burning_number) ? "" : string(burning_number)
    bound = ismissing(objective_bound_value) ? "" :
        string(round(objective_bound_value; digits = 6))
    gap = ismissing(gap_value) ? "" : string(round(gap_value; digits = 8))
    open(path, "a") do io
        println(io, join((
            graph_name,
            vertices,
            edges,
            greedy_method,
            greedy_bn,
            bn,
            bound,
            gap,
            round(time_greedy; digits = 6),
            round(time_model; digits = 6),
            round(time_solver; digits = 6),
            round(time_total; digits = 6),
            length(sources),
            status,
        ), ","))
    end
end

"""
    write_csv_result(path, graph_name; ...)

Grava um resultado completo em um arquivo temporario e depois o coloca no
caminho definitivo. Isso permite salvar a solucao da heuristica antes de
construir o PLI e substitui-la com seguranca quando o solver terminar.
"""
function write_csv_result(path::AbstractString, graph_name::AbstractString;
                          vertices::Int, edges::Int,
                          greedy_method::AbstractString, greedy_bn,
                          burning_number, objective_bound_value, gap_value,
                          time_greedy::Real, time_model::Real,
                          time_solver::Real, time_total::Real, sources, status)
    temporary_path = path * ".tmp"
    isfile(temporary_path) && rm(temporary_path)

    try
        init_csv(temporary_path)
        append_csv_row(
            temporary_path,
            graph_name;
            vertices = vertices,
            edges = edges,
            greedy_method = greedy_method,
            greedy_bn = greedy_bn,
            burning_number = burning_number,
            objective_bound_value = objective_bound_value,
            gap_value = gap_value,
            time_greedy = time_greedy,
            time_model = time_model,
            time_solver = time_solver,
            time_total = time_total,
            sources = sources,
            status = status,
        )
        mv(temporary_path, path; force = true)
    finally
        isfile(temporary_path) && rm(temporary_path)
    end
end

"Executa uma consulta ao solver e devolve `missing` se ela não estiver disponível."
function solver_value_or_missing(query, model)
    try
        value = query(model)
        return isfinite(value) ? value : missing
    catch
        return missing
    end
end

# ---------------------------------------------------------------------
# Resolução de uma instância
# ---------------------------------------------------------------------

"""
    resolve_instancia(name, g; time_limit, csv_path, warm_start, shrink_horizon)

Resolve o PLI para `g`, opcionalmente com warm start guloso e horizonte
reduzido, e grava uma linha no CSV. Se o modelo falhar ou o solver não devolver
uma solução, retorna e grava a solução válida produzida pela heurística.

Os tempos do guloso, da construção do modelo, do solver e o total são
registrados separadamente no CSV.
"""
function resolve_instancia(name::AbstractString, g::AbstractGraph;
                           csv_path::AbstractString,
                           time_limit::Real = TIME_LIMIT,
                           warm_start::Bool = WARM_START,
                           shrink_horizon::Bool = SHRINK_HORIZON,
                           silent::Bool = SILENT,
                           threads::Int = SOLVER_THREADS,
                           max_binary_variables::Int = MAX_BINARY_VARIABLES)
    start_total = time()
    n = nv(g)
    T = n  # limitante superior simples

    # O guloso sempre roda: seu valor é gravado no CSV mesmo com as flags desligadas.
    used_greedy = "grau_relativo_gatilho"
    start_greedy = time()
    z_g, src_g, states_g = greedy_2burning(g)
    time_greedy = time() - start_greedy
    println("[", name, "] guloso ", used_greedy, ": z = ", z_g,
            " | tempo = ", round(time_greedy; digits = 3), " s")
    shrink_horizon && (T = z_g)

    # x tem n(T+1) variáveis, s tem nT e f tem T.
    binary_variables = n * (T + 1) + n * T + T
    println("[", name, "] modelo: n = ", n, " | m = ", ne(g),
            " | T = ", T, " | binárias ~= ", binary_variables)

    # Resultado de seguranca. Ele e gravado antes da construcao do modelo,
    # pois essa etapa nao faz parte do limite de tempo do solver e pode falhar
    # em instancias muito grandes. Se o processo for interrompido, o CSV ainda
    # contem uma solucao valida produzida pela heuristica.
    write_csv_result(
        csv_path,
        name;
        vertices = n,
        edges = ne(g),
        greedy_method = used_greedy,
        greedy_bn = z_g,
        burning_number = z_g,
        objective_bound_value = missing,
        gap_value = missing,
        time_greedy = time_greedy,
        time_model = 0.0,
        time_solver = 0.0,
        time_total = time() - start_total,
        sources = src_g,
        status = "HEURISTIC_CHECKPOINT",
    )

    if max_binary_variables > 0 && binary_variables > max_binary_variables
        time_total = time() - start_total
        status = "HEURISTIC_FALLBACK_MODEL_TOO_LARGE"
        write_csv_result(
            csv_path,
            name;
            vertices = n,
            edges = ne(g),
            greedy_method = used_greedy,
            greedy_bn = z_g,
            burning_number = z_g,
            objective_bound_value = missing,
            gap_value = missing,
            time_greedy = time_greedy,
            time_model = 0.0,
            time_solver = 0.0,
            time_total = time_total,
            sources = src_g,
            status = status,
        )
        println("[", name, "] pulado: modelo acima do limite de ",
                max_binary_variables, " variáveis binárias")
        return z_g, src_g
    end

    start_model = time()
    local model
    local time_model
    try
        model = build_2burning_model(
            g,
            T;
            silent = silent,
            time_limit = time_limit,
            threads = threads,
        )
        warm_start && set_warm_start!(model, g, T, z_g, src_g, states_g)
        time_model = time() - start_model
        optimize!(model)
    catch e
        time_model = time() - start_model
        time_total = time() - start_total
        write_csv_result(
            csv_path,
            name;
            vertices = n,
            edges = ne(g),
            greedy_method = used_greedy,
            greedy_bn = z_g,
            burning_number = z_g,
            objective_bound_value = missing,
            gap_value = missing,
            time_greedy = time_greedy,
            time_model = time_model,
            time_solver = 0.0,
            time_total = time_total,
            sources = src_g,
            status = "HEURISTIC_FALLBACK_ERROR",
        )
        @error "Falha ao construir ou resolver o PLI de $name; foi mantida a solucao heuristica" exception = (e, catch_backtrace())
        return z_g, src_g
    end

    solver_status = termination_status(model)
    time_solver = solve_time(model)

    b2, sources = z_g, src_g
    status = solver_status
    if has_values(model)
        b2, sources = extract_solution(model, g)
    else
        status = string(solver_status, "_HEURISTIC_FALLBACK")
    end
    bound = solver_value_or_missing(objective_bound, model)
    gap = solver_value_or_missing(relative_gap, model)
    time_total = time() - start_total

    write_csv_result(
        csv_path,
        name;
        vertices = n,
        edges = ne(g),
        greedy_method = used_greedy,
        greedy_bn = z_g,
        burning_number = b2,
        objective_bound_value = bound,
        gap_value = gap,
        time_greedy = time_greedy,
        time_model = time_model,
        time_solver = time_solver,
        time_total = time_total,
        sources = sources,
        status = status,
    )
    println("[", name, "] status = ", status, " | b2 = ", b2,
            " | limite = ", bound, " | gap = ", gap,
            " | solver = ", round(time_solver; digits = 2), " s",
            " | total = ", round(time_total; digits = 2), " s")

    return b2, sources
end


"""
Percorre `input_path` recursivamente (ou processa um único arquivo), resolve
cada `.txt` e grava um CSV por grafo em `output_dir`. Pula grafos já processados,
a menos que `overwrite` seja verdadeiro.
"""
function processa_entrada(input_path::AbstractString, output_dir::AbstractString;
                          time_limit::Real = TIME_LIMIT,
                          warm_start::Bool = WARM_START,
                          shrink_horizon::Bool = SHRINK_HORIZON,
                          silent::Bool = SILENT,
                          threads::Int = SOLVER_THREADS,
                          max_binary_variables::Int = MAX_BINARY_VARIABLES,
                          overwrite::Bool = false)
    arquivos = String[]
    input_dir = isdir(input_path) ? input_path : dirname(input_path)

    if isfile(input_path)
        endswith(lowercase(input_path), ".txt") ||
            error("O arquivo de entrada precisa ter extensão .txt")
        push!(arquivos, input_path)
    else
        for (root, _, files) in walkdir(input_path)
            for f in files
                endswith(lowercase(f), ".txt") &&
                    push!(arquivos, joinpath(root, f))
            end
        end
    end

    # Instâncias menores primeiro: assim uma instância grande não impede que
    # o lote produza os resultados exatos mais fáceis.
    sort!(arquivos; by = path -> (filesize(path), path))

    total = length(arquivos)
    total > 0 || error("Nenhum arquivo .txt encontrado em $input_path")
    for (i, filepath) in enumerate(arquivos)
        rel_dir = relpath(dirname(filepath), input_dir)
        out_sub = rel_dir == "." ? output_dir : joinpath(output_dir, rel_dir)
        name    = splitext(basename(filepath))[1]
        out_csv = joinpath(out_sub, name * ".csv")

        println("[", i, "/", total, "] ", filepath)

        # Um checkpoint heurístico indica que a execução anterior foi
        # interrompida durante a construção ou resolução do modelo. Ele deve
        # ser reavaliado para receber um estado final, inclusive quando o novo
        # limite preventivo classificar o modelo como grande demais.
        resultado_finalizado = false
        if isfile(out_csv) && countlines(out_csv) >= 2
            ultima_linha = last(readlines(out_csv))
            status = strip(last(split(ultima_linha, ',')))
            resultado_finalizado = status != "HEURISTIC_CHECKPOINT"
        end

        if !overwrite && resultado_finalizado
            println("  Já processado, pulando -> ", out_csv)
            continue
        end

        mkpath(out_sub)
        try
            g = read_simple_graph(filepath)
            isfile(out_csv) && rm(out_csv)
            init_csv(out_csv)
            resolve_instancia(
                name,
                g;
                csv_path = out_csv,
                time_limit = time_limit,
                warm_start = warm_start,
                shrink_horizon = shrink_horizon,
                silent = silent,
                threads = threads,
                max_binary_variables = max_binary_variables,
            )
        catch e
            @error "Falha em $filepath" exception = (e, catch_backtrace())
            # remove CSV incompleto (só cabeçalho) para ser refeito depois
            isfile(out_csv) && countlines(out_csv) < 2 && rm(out_csv)
        end

        GC.gc()  # libera memória entre instâncias grandes
    end
end

# ---------------------------------------------------------------------
# Função principal
# ---------------------------------------------------------------------
function usage()
    println(stderr, """
Uso:
  julia --project=. burning2_pli_warmstart.jl <entrada> <saida> [opções]

Entrada pode ser um arquivo .txt ou uma pasta. Opções:
  --tempo SEGUNDOS       limite do HiGHS por grafo (padrão: 300)
  --threads N            threads do HiGHS (padrão: 1)
  --max-binarias N       pula modelos acima do limite; zero desativa
  --sem-warm-start       não fornece a solução inicial ao HiGHS
  --sem-reduzir-horizonte usa T=n em vez do valor guloso
  --log-solver           mostra o log do HiGHS
  --sobrescrever         refaz resultados já existentes
""")
end

function parse_arguments(args)
    length(args) >= 2 || (usage(); exit(1))

    input_path = abspath(args[1])
    output_dir = abspath(args[2])
    time_limit = TIME_LIMIT
    threads = SOLVER_THREADS
    max_binary_variables = MAX_BINARY_VARIABLES
    warm_start = WARM_START
    shrink_horizon = SHRINK_HORIZON
    silent = SILENT
    overwrite = false

    i = 3
    while i <= length(args)
        option = args[i]
        if option in ("--tempo", "--threads", "--max-binarias")
            i < length(args) || error("Falta o valor de $option")
            value = args[i + 1]
            if option == "--tempo"
                time_limit = parse(Float64, value)
            elseif option == "--threads"
                threads = parse(Int, value)
            elseif option == "--max-binarias"
                max_binary_variables = parse(Int, value)
            end
            i += 2
            continue
        elseif option == "--sem-warm-start"
            warm_start = false
        elseif option == "--sem-reduzir-horizonte"
            shrink_horizon = false
        elseif option == "--log-solver"
            silent = false
        elseif option == "--sobrescrever"
            overwrite = true
        elseif option in ("-h", "--help")
            usage()
            exit(0)
        else
            error("Opção desconhecida: $option")
        end
        i += 1
    end

    time_limit > 0 || error("--tempo precisa ser maior que zero")
    threads > 0 || error("--threads precisa ser maior que zero")
    max_binary_variables >= 0 || error("--max-binarias não pode ser negativo")
    ispath(input_path) || error("Entrada não encontrada: $input_path")

    return (;
        input_path,
        output_dir,
        time_limit,
        threads,
        max_binary_variables,
        warm_start,
        shrink_horizon,
        silent,
        overwrite,
    )
end

function main()
    config = parse_arguments(ARGS)
    processa_entrada(
        config.input_path,
        config.output_dir;
        time_limit = config.time_limit,
        threads = config.threads,
        warm_start = config.warm_start,
        shrink_horizon = config.shrink_horizon,
        silent = config.silent,
        max_binary_variables = config.max_binary_variables,
        overwrite = config.overwrite,
    )
end

# ---------------------------------------------------------------------
# Chamada
# ---------------------------------------------------------------------
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
