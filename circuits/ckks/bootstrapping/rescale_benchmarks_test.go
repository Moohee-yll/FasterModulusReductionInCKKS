package bootstrapping

import (
	"fmt"
	"testing"

	"github.com/stretchr/testify/require"

	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	"github.com/tuneinsight/lattigo/v6/utils"
	"github.com/tuneinsight/lattigo/v6/utils/sampling"
)

func BenchmarkRescale(b *testing.B) {
	schemeParamsLit := test_rescale_set1.SchemeParams
	btpParamsLit := ParametersLiteral{}

	if *flagLongTest {
		schemeParamsLit.LogN = 16
	}

	params, err := ckks.NewParametersFromLiteral(schemeParamsLit)
	require.Nil(b, err)

	btpParamsLit.LogN = utils.Pointy(params.LogN())

	btpParams, err := NewParametersFromLiteral(params, btpParamsLit)
	require.Nil(b, err)

	// Insecure params for fast testing only
	if !*flagLongTest {
		btpParams.SlotsToCoeffsParameters.LogSlots = btpParams.BootstrappingParameters.LogN() - 1
		btpParams.CoeffsToSlotsParameters.LogSlots = btpParams.BootstrappingParameters.LogN() - 1

		// Corrects the message ratio to take into account the smaller number of slots and keep the same precision
		btpParams.Mod1ParametersLiteral.LogMessageRatio += 16 - params.LogN()
	}
	b.Logf("ParamsN2: LogN=%d/LogSlots=%d/LogQP=%f", params.LogN(), params.LogMaxSlots(), params.LogQP())

	// sk := rlwe.NewKeyGenerator(btpParams.BootstrappingParameters).GenSecretKeyNew()

	b.Log("Generating Bootstrapping Keys")
	// btpKeys, _, err := btpParams.GenEvaluationKeys(sk)
	// require.NoError(t, err)

	// evaluator, err := NewEvaluator(btpParams, btpKeys)
	// require.NoError(t, err)

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

	// ctQ0, err := enc.EncryptNew(plaintext)

	require.NoError(b, err)

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
	// warmup := 100

	b.Run("ModDOwn", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			be.ModDownQPtoQ(polyQ.Level(), polyP.Level(), polyQ, polyP, polyQ2)
		}
	})

	fmt.Println("qk: ", ringQP.RingP.ModuliChainLength())

	b.Run(fmt.Sprint("Rescale"), func(b *testing.B) {
		// for i := 0; i < b.N; i++ {
		for i := 0; i < 200; i++ {
			for idx, _ := range ringQP.RingP.ModuliChain() {

				ringQP_singlering.AtLevel(rQPLvl-idx).DivRoundByLastModulusAux(polyQP, polyQP2, 0)
			}
		}

	})

	rpLvl := polyP.Level()
	fmt.Println("rplvl: ", rpLvl)

	b.Run(fmt.Sprint("Rescale Optimized"), func(b *testing.B) {
		lvl := rpLvl + 1
		moduli := make([]uint64, lvl)
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
				// scalarMonts[j] = v
			}
		}

		for i := 0; i < b.N; i++ {
			for lvl := 0; lvl < rpLvl; lvl++ {
				ringQP.RingP.AtLevel(rpLvl-lvl).DivFloorByLastModulus(polyP, polyP)
			}

			ringQP_singlering.DivRoundByLastModulusKernelOpt(polyQP, rpLvl+1, moduli, medconsts, scalarMonts)
			// for idx, _ := range ringQP.RingP.ModuliChain() {

			// 	ringQP_singlering.AtLevel(rQPLvl-idx).DivRoundByLastModulusAux(polyQP, polyQP2, 0)
			// 	// polyQP2.Resize(rQPLvl - idx - 1)
			// }
		}

	})
}
