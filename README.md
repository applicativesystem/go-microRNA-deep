# go-microRNAs-deep

- a microRNA go package for implementing the deep learning approaches.
- takes the predicted microRNAs, fasta and extracts the fasta, upstream, downstream, tokenization and neural network. 
- prepares the microRNA predictions from all microRNA target predictions tools into a single package and gives you the structured data with tokenization for direct input as weights into the neural networks. 
- I implemented this package in both GO and RUST for a special issue in Springer Invitation. 

```
gauavsablok@gauravsablok ~/Desktop/go/go-microRNA-deep ±main⚡ » \
go run main.go
This analyzes the microRNA prediction and makes them ready for the deep learning approaches

Usage:
  analyzePred [command]

Available Commands:
  completion      Generate the autocompletion script for the specified shell
  help            Help about any command
  psRNAanalyzer
  psRNAmapanalyze
  tapiranalyzer
  tarHunter
  targetFinder

Flags:
  -h, --help   help for analyzePred

Use "analyzePred [command] --help" for more information about a command.
gauavsablok@gauravsablok ~/Desktop/go/go-microRNA-deep ±main⚡ » \
go run main.go psRNAanalyzer -h
Analyzes and prepares the psRNA target predictions for the deep learning

Usage:
  analyzePred psRNAanalyzer [flags]

Flags:
      --expectation value float   expectation value (default 0.5)
  -f, --fastapred string          fasta predict (default "fasta file for the predictions")
  -h, --help                      help for psRNAanalyzer
  -p, --psRNAPred string          psRNA predictions (default "psRNA microRNA predictions")
gauavsablok@gauravsablok ~/Desktop/go/go-microRNA-deep ±main⚡ » \
go run main.go psRNAmapanalyze -h
Analyze the psRNA map alignment of the reads to the genome

Usage:
  analyzePred psRNAmapanalyze [flags]

Flags:
  -f, --fastapred string   fasta predict (default "fasta file for the predictions")
  -h, --help               help for psRNAmapanalyze
  -P, --psRNAfile string   map reads to the genome file (default "RNA mapping file")
gauavsablok@gauravsablok ~/Desktop/go/go-microRNA-deep ±main⚡ » \
go run main.go tapiranalyzer -h                                                                    
Analyzes tapir target predictions for the deep learning

Usage:
  analyzePred tapiranalyzer [flags]

Flags:
  -f, --fastapred string   fasta predict (default "fasta file for the predictions")
  -h, --help               help for tapiranalyzer
  -p, --psRNAPred string   psRNA predictions (default "psRNA microRNA predictions")
gauavsablok@gauravsablok ~/Desktop/go/go-microRNA-deep ±main⚡ » \
go run main.go tarHunter -h
analyze the tarHunter results for the miRNA predictions

Usage:
  analyzePred tarHunter [flags]

Flags:
  -D, --downstream int         downstream of the miRNA predictions (default 10)
  -f, --fastapred string       fasta predict (default "fasta file for the predictions")
  -h, --help                   help for tarHunter
  -T, --tarhunterfile string   tarhunter analysis (default "tarHunter predictions")
  -U, --upstream int           upstream of the miRNA predictions (default 10)
gauavsablok@gauravsablok ~/Desktop/go/go-microRNA-deep ±main⚡ » \
go run main.go targetFinder -h
analyzes the targetFinder results for the miRNA predictions

Usage:
  analyzePred targetFinder [flags]

Flags:
  -D, --downstream int            downstream of the miRNA predictions (default 10)
  -f, --fastapred string          fasta predict (default "fasta file for the predictions")
  -h, --help                      help for targetFinder
  -T, --targetFinderfile string   targetFinder analysis (default "targetFinder predictions")
  -U, --upstream int              upstream of the miRNA predictions (default 10)
``

Gaurav Sablok
