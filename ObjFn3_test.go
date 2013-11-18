package ser

// Test of obj fns aganist Kostopoulos values
import (
	"fmt"
	"os"
	"testing"
)

func TestObjFns2(t *testing.T) {
	file, err := os.Open("/home/pac/live/ext/src/code.google.com/p/ser/dat/artif/matrix/randR2.csv") // For read access.
	if err != nil {
		panic("file does not exist")
	}
	a := ReadCsvMatrix64(file)

	n, m := a.Dims()
	p1 := NewIntVector(n)
	p1.Order()
	p1.Print()
	p2 := NewIntVector(m)
	p2.Order()
	p2.Print()
	//rectangular matrix
	fmt.Println("MEffGain: ", MEffGain(a, p1, p2))
	fmt.Println("MooreStressLoss: ", MooreStressLoss(a, p1, p2))
	fmt.Println("VonNeumannStressLoss: ", VonNeumannStressLoss(a, p1, p2))
	fmt.Println("BertinGain: ", BertinLoss(a, p1, p2))
}

/*
moe:  -4.1300
mstr: 1.4436e+03
nstr: 753.7555
bcc:  15.1998

MEffGain:  -4.130017559673509
MooreStressLoss:  1443.6377804913923
VonNeumannStressLoss:  753.7555028079539
BertinLoss:  -2.8526416837375566
*/
