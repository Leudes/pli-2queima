# ILP_2burning

Programa Linear Inteiro (PLI) para o **número de 2-queima** (*2-burning number*, $b_2(G)$) de grafos, implementado em Julia com [JuMP](https://jump.dev/) e [HiGHS](https://highs.dev/). O modelo recebe como solução inicial (*warm start*) o resultado de um algoritmo guloso, que também fornece um limitante superior para reduzir o horizonte de tempo.

## O processo de 2-queima

O processo de 2-queima é uma variante do processo de queima de grafos (*graph burning*) em que um vértice precisa de **dois** vizinhos queimados para queimar por propagação. O processo ocorre em rodadas $t = 1, 2, \dots$:

1. **Propagação:** todo vértice não queimado com pelo menos 2 vizinhos queimados (no estado da rodada $t-1$) passa a ficar queimado na rodada $t$.
2. **Nova fonte:** se ainda restarem vértices não queimados, escolhe-se um deles como nova fonte, que é queimado.

O processo termina na primeira rodada em que todos os vértices estão queimados. O número de 2-queima $b_2(G)$ é o menor número de rodadas necessário, minimizado sobre todas as escolhas de fontes.

## Modelo PLI

Variáveis (com $T$ o horizonte de tempo):

| Variável | Significado |
|---|---|
| `x[v,t]` | 1 se o vértice `v` está queimado até o instante `t` ($t = 0,\dots,T$) |
| `s[v,t]` | 1 se `v` é a fonte escolhida no instante `t` ($t = 1,\dots,T$) |
| `f[t]` | 1 se o processo termina exatamente no instante `t` |
| `z` | instante final (valor objetivo) |

O objetivo é minimizar `z`, sujeito às restrições R1–R11 (monotonicidade de `x`, no máximo uma fonte por instante, cada vértice é fonte no máximo uma vez, regras de propagação, instante final único e todos queimados no instante final). A descrição detalhada de cada restrição está comentada no código-fonte, em `build_2burning_model`.

### Algoritmo guloso e warm start

A solução inicial combina o grau relativo com a quantidade de gatilhos criados, usando exatamente a mesma regra da solução inicial oficial do GVNS. O valor obtido, $z_g$, é um limitante superior válido para $b_2(G)$ e é usado de duas formas:

- **`WARM_START`**: a solução gulosa completa (`x`, `s`, `f`, `z`) é fornecida ao HiGHS como solução inicial (MIP start).
- **`SHRINK_HORIZON`**: o horizonte passa a ser $T = z_g$ em vez de $T = n$, o que reduz bastante o número de variáveis e restrições.

## Requisitos

- [Julia](https://julialang.org/downloads/) 1.10 ou superior
- Pacotes: `JuMP`, `Graphs`, `HiGHS`

O HiGHS é instalado automaticamente pelo gerenciador de pacotes do Julia e não exige licença nem instalação externa.

Instalação dos pacotes Julia no ambiente do projeto:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

## Uso

```bash
julia --project=. burning2_pli_warmstart.jl <entrada> <dir_saida> [opções]
```

Exemplo (é o que o script `run.sh` executa):

```bash
./run.sh
```

Para executar diretamente a base final deste TCC, que fica no repositório
vizinho `BRKGA`:

```bash
./run_base_final.sh
```

As principais opções são `--tempo SEGUNDOS`, `--threads N`,
`--log-solver` e `--sobrescrever`. A entrada pode ser uma pasta ou apenas um
arquivo `.txt`.
O script da base usa `--max-binarias 0`, portanto tenta construir o modelo de
todas as instâncias. Se uma instância falhar antes de o solver produzir uma
solução, o CSV conserva como solução viável de segurança o resultado da mesma
heurística inicial usada pelo GVNS.

O programa percorre `<dir_entrada>` recursivamente, resolve cada arquivo `.txt` e grava **um CSV por grafo** em `<dir_saida>`, espelhando a estrutura de subpastas da entrada:

```
input/                      output/
├── hb/                     ├── hb/
│   └── can__96.txt         │   └── can__96.csv
└── paths/                  └── paths/
    ├── p5.txt                  ├── p5.csv
    ├── p6.txt                  ├── p6.csv
    └── p15.txt                 └── p15.csv
```

Todas as instâncias são resolvidas em uma única sessão Julia, o que evita pagar o tempo de carga e compilação do JuMP/HiGHS a cada grafo. Grafos que já possuem CSV completo (cabeçalho e pelo menos uma linha) são **pulados**, então é possível interromper e retomar a execução. Se a resolução de um grafo falhar, o erro é registrado e o lote continua com os demais.

## Formato de entrada

Arquivos `.txt` com uma aresta por linha, dada por dois inteiros separados por espaço:

```
0 1
1 2
2 3
3 4
```

- Os rótulos dos vértices podem ser quaisquer inteiros; internamente são normalizados para índices consecutivos.
- Laços e arestas repetidas são descartados (o grafo é simples e não direcionado).
- Linhas vazias e linhas iniciadas por `#` ou `%` são ignoradas, assim como colunas além da segunda.

## Formato de saída

Cada CSV tem o cabeçalho e uma linha:

| Coluna | Descrição |
|---|---|
| `graph` | nome do grafo (nome do arquivo sem extensão) |
| `vertices`, `edges` | tamanho da instância lida |
| `greedy_method` | guloso usado para construir a solução inicial |
| `greedy_burning_number` | valor obtido pelo algoritmo guloso |
| `pli_burning_number` | melhor valor encontrado pelo PLI (vazio se nenhuma solução viável foi encontrada) |
| `objective_bound` | melhor limitante inferior provado pelo HiGHS |
| `gap` | distância relativa entre a melhor solução e o limitante |
| `time_greedy` | tempo da solução inicial |
| `time_model` | tempo para construir o modelo e inserir o *warm start* |
| `time_solver` | tempo informado pelo HiGHS |
| `time_total` | tempo total da instância, sem a leitura do arquivo |
| `#sources` | número de fontes na solução do PLI |
| `status` | status de terminação do solver (`OPTIMAL`, `TIME_LIMIT`, ...) |

Exemplo (`output/paths/p15.csv`):

```csv
graph,vertices,edges,greedy_method,greedy_burning_number,pli_burning_number,objective_bound,gap,time_greedy,time_model,time_solver,time_total,#sources,status
p15,15,14,grau_relativo_gatilho,9,9,9.0,0.0,0.001,0.5,5.522,6.023,9,OPTIMAL
```

> Quando o `status` é `TIME_LIMIT`, o valor em `pli_burning_number` é a melhor solução encontrada dentro do limite de tempo e **não** é necessariamente ótimo. No exemplo incluído, `can__96` termina com `TIME_LIMIT`, com guloso 10 e PLI 8.

## Configuração

Os parâmetros ficam no início de `burning2_pli_warmstart.jl`:

| Constante | Padrão | Descrição |
|---|---|---|
| `TIME_LIMIT` | `300.0` | limite de tempo do solver por grafo, em segundos |
| `WARM_START` | `true` | fornece a solução gulosa como solução inicial ao HiGHS |
| `SHRINK_HORIZON` | `true` | usa $T = z_g$ em vez de $T = n$ |
| `SILENT` | `true` | `false` exibe o log do HiGHS (útil para conferir o MIP start) |
| `SOLVER_THREADS` | `1` | quantidade padrão de threads do solver |
| `MAX_BINARY_VARIABLES` | `0` | limite preventivo; zero desativa |

Para rodar vários processos em paralelo, convém limitar as threads do HiGHS para que as instâncias não disputem CPU e os tempos fiquem comparáveis.

## Estrutura do repositório

```
ILP_2burning/
├── burning2_pli_warmstart.jl   # modelo PLI, guloso, leitura de grafos e execução em lote
├── Project.toml                # dependências Julia reproduzíveis
├── Manifest.toml               # versões exatas das dependências instaladas
├── run.sh                      # exemplo de execução
├── run_base_final.sh           # executa a base final do TCC
├── input/                      # grafos de exemplo (listas de arestas)
├── output/                     # CSVs gerados
├── LICENSE
└── README.md
```

## Licença

Distribuído sob a licença [GNU GPL v3](LICENSE).

## Como citar

Se você usar este código em um trabalho acadêmico, considere citar o repositório (adicione aqui a referência ou o DOI quando houver).
