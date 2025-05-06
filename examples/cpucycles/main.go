package main

// #include "tsc.c"
import "C"

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v6/circuits/ckks/bootstrapping"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	"github.com/tuneinsight/lattigo/v6/utils"
	"github.com/tuneinsight/lattigo/v6/utils/sampling"
)

func main() {
	schemeParamsLit := ckks.ParametersLiteral{
		LogN:            16,
		LogQ:            []int{60, 40, 40, 40, 40, 40, 40, 40},
		LogP:            []int{61, 61, 61},
		Xs:              ring.Ternary{H: 192},
		LogDefaultScale: 40,
	}

	result_moddown := [][]C.uint64_t{}
	result_rescale := [][]C.uint64_t{}
	result_opt_rescale := [][]C.uint64_t{}
	for numP := 4; numP <= 7; numP++ {
		result_ki_moddown := []C.uint64_t{}
		result_ki_rescale := []C.uint64_t{}
		result_ki_opt_rescale := []C.uint64_t{}
		schemeParamsLit.LogP = append(schemeParamsLit.LogP, 61)
		schemeParamsLit.LogQ = []int{60} //{60, 40, 40, 40, 40, 40, 40}
		for numQ := 1; numQ <= 31; numQ++ {
			fmt.Println("numQ,", numQ, ", numP,", numP)
			res1, res2, res3 := count_cycles(schemeParamsLit)
			result_ki_moddown = append(result_ki_moddown, res1)
			result_ki_rescale = append(result_ki_rescale, res2)
			result_ki_opt_rescale = append(result_ki_opt_rescale, res3)
			schemeParamsLit.LogQ = append(schemeParamsLit.LogQ, 40)
		}
		result_moddown = append(result_moddown, result_ki_moddown)
		result_rescale = append(result_rescale, result_ki_rescale)
		result_opt_rescale = append(result_opt_rescale, result_ki_opt_rescale)
	}
	fmt.Println("result_moddown", result_moddown)
	fmt.Println("result_rescale", result_rescale)
	fmt.Println("result_opt_rescale", result_opt_rescale)
}

func count_cycles(schemeParamsLit ckks.ParametersLiteral) (C.uint64_t, C.uint64_t, C.uint64_t) {

	test_round := 100

	btpParamsLit := bootstrapping.ParametersLiteral{}

	flagLongTest := true

	if flagLongTest {
		schemeParamsLit.LogN = 16
	}

	params, err := ckks.NewParametersFromLiteral(schemeParamsLit)
	if err != nil {
		panic(err)
	}

	btpParamsLit.LogN = utils.Pointy(params.LogN())

	btpParams, err := bootstrapping.NewParametersFromLiteral(params, btpParamsLit)

	// Insecure params for fast testing only
	if !flagLongTest {
		btpParams.SlotsToCoeffsParameters.LogSlots = btpParams.BootstrappingParameters.LogN() - 1
		btpParams.CoeffsToSlotsParameters.LogSlots = btpParams.BootstrappingParameters.LogN() - 1

		// Corrects the message ratio to take into account the smaller number of slots and keep the same precision
		btpParams.Mod1ParametersLiteral.LogMessageRatio += 16 - params.LogN()
	}

	be := ring.NewBasisExtender(btpParams.ResidualParameters.RingQ(), btpParams.ResidualParameters.RingP())
	ringQP := btpParams.ResidualParameters.RingQP()

	values := make([]complex128, params.MaxSlots())
	for i := range values {
		values[i] = sampling.RandComplex128(-1, 1)
	}

	values[0] = complex(0.9238795325112867, 0.3826834323650898)
	values[1] = complex(0.9238795325112867, 0.3826834323650898)
	if len(values) > 2 {
		values[2] = complex(0.9238795325112867, 0.3826834323650898)
		values[3] = complex(0.9238795325112867, 0.3826834323650898)
	}

	ecd := ckks.NewEncoder(params)
	// enc := rlwe.NewEncryptor(params, sk)
	// dec := rlwe.NewDecryptor(params, sk)

	plaintext := ckks.NewPlaintext(params, 0)
	ecd.Encode(values, plaintext)

	prng, _ := sampling.NewKeyedPRNG([]byte{'a', 'b', 'c'})

	samplerQ := ring.NewUniformSampler(prng, ringQP.RingQ)
	samplerP := ring.NewUniformSampler(prng, ringQP.RingP)

	polyP := samplerP.ReadNew()
	polyQ := samplerQ.ReadNew()
	// levelP := polyP.Level()
	fmt.Println("polyQ: ", polyQ.Level(), "polyP: ", polyP.Level())

	polyQ2 := ringQP.RingQ.NewPoly()
	moduli := append(ringQP.RingQ.ModuliChain(), ringQP.RingP.ModuliChain()...)
	ringQP_singlering, _ := ring.NewRing(ringQP.N(), moduli)
	samplerQP := ring.NewUniformSampler(prng, ringQP_singlering)
	polyQP2 := ringQP_singlering.NewPoly()
	// polyQP3 := ringQP_singlering.NewPoly()
	polyQP := samplerQP.ReadNew()
	rQPLvl := ringQP_singlering.Level()
	// num_run := 100
	warmup := 100

	// warming up
	for i := 0; i < warmup; i++ {
		be.ModDownQPtoQ(polyQ.Level(), polyP.Level(), polyQ, polyP, polyQ2)
	}

	moddown_start := C.tsc()
	for i := 0; i < test_round; i++ {
		be.ModDownQPtoQ(polyQ.Level(), polyP.Level(), polyQ, polyP, polyQ2)
	}
	moddown_end := C.tsc()
	moddown_cycles := (moddown_end - moddown_start) / C.uint64_t(test_round)
	fmt.Println("ModDown cycles:", (moddown_end-moddown_start)/C.uint64_t(test_round))

	for i := 0; i < warmup; i++ {
		for idx, _ := range ringQP.RingP.ModuliChain() {
			ringQP_singlering.AtLevel(rQPLvl-idx).DivRoundByLastModulusAux(polyQP, polyQP2, 0)
		}
	}

	rescale_start := C.tsc()
	for i := 0; i < test_round; i++ {
		for idx, _ := range ringQP.RingP.ModuliChain() {
			ringQP_singlering.AtLevel(rQPLvl-idx).DivRoundByLastModulusAux(polyQP, polyQP2, 0)
		}
	}
	rescale_end := C.tsc()
	rescale_cycles := (rescale_end - rescale_start) / C.uint64_t(test_round)
	fmt.Println("Rescale cycles:", (rescale_end-rescale_start)/C.uint64_t(test_round))

	rpLvl := polyP.Level()
	lvl := rpLvl + 1
	moduli = make([]uint64, lvl)
	medconsts := make([]uint64, lvl)
	level := ringQP_singlering.Level()
	for j, s := range ringQP_singlering.SubRings[level+1-lvl:] {
		moduli[j] = s.Modulus
		medconsts[j] = s.MRedConstant
	}
	scalarMonts := make([][]uint64, level+1-lvl)

	for i, _ := range scalarMonts {
		scalarMonts[i] = make([]uint64, lvl)
		for j, _ := range ringQP_singlering.SubRings[level+1-lvl:] {
			scalarMonts[i][j] = ringQP_singlering.RescaleConstants[level-1-j][i]
		}
	}

	for i := 0; i < warmup; i++ {
		for lvl := 0; lvl < rpLvl; lvl++ {
			ringQP.RingP.AtLevel(rpLvl-lvl).DivFloorByLastModulus(polyP, polyP)
		}

		ringQP_singlering.DivRoundByLastModulusKernelOpt(polyQP, rpLvl+1, moduli, medconsts, scalarMonts)
	}
	res_start := C.tsc()
	for i := 0; i < test_round; i++ {
		for lvl := 0; lvl < rpLvl; lvl++ {
			ringQP.RingP.AtLevel(rpLvl-lvl).DivFloorByLastModulus(polyP, polyP)
		}

		ringQP_singlering.DivRoundByLastModulusKernelOpt(polyQP, rpLvl+1, moduli, medconsts, scalarMonts)
	}
	res_end := C.tsc()

	fmt.Println("rescale opt cycles:", (res_end-res_start)/C.uint64_t(test_round))
	rescale_opt_cycles := (res_end - res_start) / C.uint64_t(test_round)
	return moddown_cycles, rescale_cycles, rescale_opt_cycles
}
