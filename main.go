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
	psRNAfile string
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

var psRNAMapCmd = &cmd.Command{
	Use:  "psRNAmapanalyze",
	Long: "Analyze the psRNA map alignment of the reads to the genome",
	Run:  psRNAMapFunc,
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
	psRNAMapCmd.Flags().
		StringVarP(&psRNAfile, "psRNAfile", "P", "RNA mapping file", "map reads to the genome file")
	psRNAMapCmd.Flags().
		StringVarP(&fastPred, "fastapred", "f", "fasta file for the predictions", "fasta predict")

	rootCmd.AddCommand(psRNACmd)
	rootCmd.AddCommand(tapirCmd)
	rootCmd.AddCommand(psRNAMapCmd)
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

func tapirFunc(cmd *cobra.Command, args []string) {
	target := []string{}
	miRNA := []string{}
	score := []int{}
	mfe := []int{}
	start := []int{}
	end := []int{}

	tapirOpen, err := os.Open(tapirPred)
	if err != nil {
		log.Fatal(err)
	}
	tapirRead := bufio.NewScanner(tapirOpen)
	for tapirRead.Scan() {
		line := tapirRead.Text()
		if strings.HasPrefix(string(line), "target") {
			target = append(target, strings.Split(string(line), " ")[1])
		}
		if strings.HasPrefix(string(line), "miRNA") {
			miRNA := strings.Split(strings.Split(string(line), ":")[3], "=")[2]
		}
		if strings.HasPrefix(string(line), "score") {
			scoreInter, _ := strconv.Atoi(strings.Split(string(line), " ")[1])
			score = append(score, scoreInter)
		}
		if strings.HasPrefix(string(line), "mfe") {
			mfeInter, _ := strconv.Atoi(strings.Split(string(line), " ")[1])
			mfe = append(mfe, mfeInter)
		}
		if strings.HasPrefix(string(line), "start") {
			startInter, _ := strconv.Atoi(strings.Split(string(line), " ")[1])
			start = append(start, startInter)
		}
		if strings.HasPrefix(string(line), "target_5'") {
			end = append(end, len(strings.Split(string(line), " ")[1]))
		}
	}

	type tapirExtract struct {
		tapid  string
		tapseq string
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

	tapSeq := []tapirExtract{}

	for i := range target {
		for j := range fastID {
			if target[i] == fastID[i].id {
				tapSeq = append(tapSeq, tapirExtract{
					tapid:  target[i],
					tapseq: fastSeq[i].seq[start[i]:end[i]],
				})
			}
		}
	}

	tapwrite, err := os.Create("tapirneural.fasta")
	if err != nil {
		log.Fatal(err)
	}
	defer tapwrite.Close()
	for i := range tapSeq {
		tapwrite.WriteString(tapSeq[i].tapid + "\t" + tapSeq[i].tapseq)
	}
}

func psRNAMapFunc(cmd *cobra.Command, args []string) {
	type psRNAfunc struct {
		id     string
		ref    string
		strand string
		start  int
		stop   int
		read   string
	}

	readMap := []psRNAfunc{}
	fOpen, err := os.Open(psRNAfile)
	if err != nil {
		log.Fatal(err)
	}
	fRead := bufio.NewScanner(fOpen)
	for fRead.Scan() {
		line := fRead.Text()
		startStore, _ := strconv.Atoi(strings.Split(string(line), " ")[3])
		endStore, _ := strconv.Atoi(strings.Split(string(line), " ")[4])
		readMap = append(readMap, psRNAfunc{
			id:     strings.Split(string(line), " ")[0],
			ref:    strings.Split(string(line), " ")[1],
			strand: strings.Split(string(line), " ")[2],
			start:  startStore,
			stop:   endStore,
			read:   strings.Split(string(line), " ")[5],
		})
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

	type readMapSeq struct {
		read   string
		start  int
		stop   int
		refseq string
	}

	readMapSeqStore := []readMapSeq{}

	for i := range readMap {
		for j := range fastID {
			if readMap[i].id == fastID[i].id {
				readMapSeqStore = append(readMapSeqStore, readMapSeq{
					read:   readMap[i].read,
					start:  readMap[i].start,
					stop:   readMap[i].stop,
					refseq: fastSeq[j].seq[readMap[i].start:readMap[i].stop],
				})
			}
		}
	}

	readMapWrite, err := os.Create("psRNAMap.fasta")
	if err != nil {
		log.Fatal(err)
	}
	for i := range readMapSeqStore {
		intStartWrite := string(readMapSeqStore[i].start)
		intEndWrite := string(readMapSeqStore[i].stop)
		readMapWrite.WriteString(
			readMapSeqStore[i].read + "\t" + readMapSeqStore[i].refseq + "\t" + intStartWrite + "\t" + intEndWrite,
		)
	}
}