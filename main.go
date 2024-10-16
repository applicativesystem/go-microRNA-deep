package main

/*

Author Gaurav Sablok
Universitat Potsdam
Date 2024-10-16

a miRNA serialization for the extraction of the microRNA predictions and also for the filtering of the microRNAs
and prepairing the microRNAs for the deep learning. It takes all forms of the micrRNA predictions and all the target
analyzer for the microRNA predictions and prepares them for the deep learning.


*/

import (
	"bufio"
	"log"
	"os"
	"strconv"
	"strings"

	"github.com/spf13/cobra"
)

func main() {
	if err := rootCmd.Execute(); err != nil {
		log.Fatal(err)
		os.Exit(1)
	}
}

var (
	psRNAPred string
	tapirPred string
	fastPred  string
	evalue    float64
)

// var fastasequence string
// var filter float32
var rootCmd = &cobra.Command{
	Use:  "analyzePred",
	Long: "This analyzes the microRNA prediction and makes them ready for the deep learning approaches",
}

var psRNACmd = &cobra.Command{
	Use:  "psRNAanalyzer",
	Long: "Analyzes and prepares the psRNA target predictions for the deep learning",
	Run:  psRNAFunc,
}

var tapirCmd = &cobra.Command{
	Use:  "tapiranalyzer",
	Long: "Analyzes tapir target predictions for the deep learning",
	Run:  tapirFunc,
}

func init() {
	psRNACmd.Flags().
		StringVarP(&psRNAPred, "psRNAPred", "p", "psRNA microRNA predictions", "psRNA predictions")
	psRNACmd.Flags().
		StringVarP(&fastPred, "fastapred", "f", "fasta file for the predictions", "fasta predict")
	psRNACmd.Flags().
		Float64VarP(&evalue, "expectation", "expectation value for filtering", 1.5, "expectation value")
	tapirCmd.Flags().
		StringVarP(&psRNAPred, "psRNAPred", "p", "psRNA microRNA predictions", "psRNA predictions")
	tapirCmd.Flags().
		StringVarP(&fastPred, "fastapred", "f", "fasta file for the predictions", "fasta predict")
}

func psRNAFunc(cmd *cobra.Command, args []string) {
	type psRNAStruct struct {
		miRNA         string
		target        string
		evalue        float64
		targetStart   int
		targetEnd     int
		miRNAAligned  string
		targetAligned string
	}

	type psRNAStructFiltered struct {
		miRNA         string
		target        string
		evalue        float64
		targetStart   int
		targetEnd     int
		miRNAAligned  string
		targetAligned string
	}

	type fastPredID struct {
		id string
	}

	type fastPredSeq struct {
		seq string
	}

	fastID := []fastPredID{}
	fastSeq := []fastPredSeq{}

	fastaOpen, err := os.Open(fastPred)
	if err != nil {
		log.Fatal(err)
	}
	fastaRead := bufio.NewScanner(fastaOpen)
	for fastaRead.Scan() {
		line := fastaRead.Text()
		if strings.HasPrefix(string(line), ">") {
			fastID = append(fastID, fastPredID{
				id: strings.ReplaceAll(string(line), ">", ""),
			})
		}
		if !strings.HasPrefix(string(line), ">") {
			fastSeq = append(fastSeq, fastPredSeq{
				seq: string(line),
			})
		}
	}

	storemiRNA := []psRNAStruct{}
	filteredmiRNA := []psRNAStructFiltered{}

	fOpen, err := os.Open(psRNAPred)
	if err != nil {
		log.Fatal(err)
	}
	fRead := bufio.NewScanner(fOpen)
	for fRead.Scan() {
		line := fRead.Text()
		if strings.HasPrefix(string(line), "#") {
			continue
		}
		if !strings.HasPrefix(string(line), "#") {
			evalueMid, _ := strconv.ParseFloat(strings.Split(string(line), " ")[2], 64)
			targetStartMid, _ := strconv.Atoi(strings.Split(string(line), " ")[6])
			targetEndMid, _ := strconv.Atoi(strings.Split(string(line), " ")[7])
			storemiRNA = append(storemiRNA, psRNAStruct{
				miRNA:         strings.Split(string(line), " ")[0],
				target:        strings.Split(string(line), " ")[1],
				evalue:        evalueMid,
				targetStart:   targetStartMid,
				targetEnd:     targetEndMid,
				miRNAAligned:  strings.Split(string(line), " ")[8],
				targetAligned: strings.Split(string(line), " ")[9],
			})
		}
	}

	for i := range storemiRNA {
		if storemiRNA[i].evalue <= evalue {
			filteredmiRNA = append(filteredmiRNA, psRNAStructFiltered{
				miRNA:         storemiRNA[i].miRNA,
				target:        storemiRNA[i].target,
				evalue:        storemiRNA[i].evalue,
				targetStart:   storemiRNA[i].targetStart,
				targetEnd:     storemiRNA[i].targetEnd,
				miRNAAligned:  storemiRNA[i].miRNAAligned,
				targetAligned: storemiRNA[i].targetAligned,
			})
		}
	}

	type extractSeq struct {
		target      string
		targetSeq   string
		targetStart int
		targetEnd   int
	}

	targetExtract := []extractSeq{}

	for i := range filteredmiRNA {
		for j := range fastID {
			if filteredmiRNA[i].target == fastID[j].id {
				targetExtract = append(targetExtract, extractSeq{
					target:      filteredmiRNA[i].target,
					targetSeq:   fastSeq[j].seq[filteredmiRNA[i].targetStart:filteredmiRNA[i].targetEnd],
					targetStart: filteredmiRNA[i].targetStart,
					targetEnd:   filteredmiRNA[i].targetEnd,
				})
			}
		}
	}

	psRNAneural, err := os.Create("psRNANeural.txt")
	if err != nil {
		log.Fatal(err)
	}
	defer psRNAneural.Close()
	for i := range targetExtract {
		psRNAneural.WriteString(
			targetExtract[i].target + "\t" + ">" + targetExtract[i].targetSeq,
		)
	}
}
