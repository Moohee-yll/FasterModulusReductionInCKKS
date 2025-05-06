package bootstrapping

import (
	"flag"
	"fmt"
	"math"
	"testing"
	"time"

	"github.com/stretchr/testify/require"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	"github.com/tuneinsight/lattigo/v6/utils"
	"github.com/tuneinsight/lattigo/v6/utils/sampling"
)

var flagLongTest = flag.Bool("long", true, "run the long test suite (all parameters + secure bootstrapping). Overrides -short and requires -timeout=0.")
var printPrecisionStats = flag.Bool("print-precision", false, "print precision stats")

var testPrec45 = ckks.ParametersLiteral{
	LogN:            10,
	LogQ:            []int{60, 40},
	LogP:            []int{61},
	LogDefaultScale: 40,
}

func TestCtSize(t *testing.T) {
	t.Run("CtSize", func(t *testing.T) {
		schemeParamsLit :=
			ckks.ParametersLiteral{
				LogN:            16,
				LogQ:            []int{60, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 60},
				LogP:            []int{61, 61, 61, 61, 61, 61},
				Xs:              ring.Ternary{H: 32768},
				LogDefaultScale: 40,
			}

		btpParamsLit := ParametersLiteral{}

		if *flagLongTest {
			schemeParamsLit.LogN = 16
		}

		params, err := ckks.NewParametersFromLiteral(schemeParamsLit)
		require.Nil(t, err)

		btpParamsLit.LogN = utils.Pointy(params.LogN())

		btpParams, err := NewParametersFromLiteral(params, btpParamsLit)
		require.Nil(t, err)
		sk := rlwe.NewKeyGenerator(btpParams.BootstrappingParameters).GenSecretKeyNew()

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
		enc := rlwe.NewEncryptor(params, sk)

		plaintext := ckks.NewPlaintext(params, 12)
		ecd.Encode(values, plaintext)

		ctQ0, err := enc.EncryptNew(plaintext)
		fmt.Println("ctQ0: ", ctQ0.Level())
		data, _ := ctQ0.MarshalBinary()

		fmt.Println("datalen: ", len(data))
	})
}

func TestRescaling(t *testing.T) {
	t.Run("RescalingLevelByLevel", func(t *testing.T) {
		schemeParamsLit := test_rescale_set1.SchemeParams
		btpParamsLit := ParametersLiteral{}

		if *flagLongTest {
			schemeParamsLit.LogN = 16
		}

		params, err := ckks.NewParametersFromLiteral(schemeParamsLit)
		require.Nil(t, err)

		btpParamsLit.LogN = utils.Pointy(params.LogN())

		btpParams, err := NewParametersFromLiteral(params, btpParamsLit)
		require.Nil(t, err)

		// Insecure params for fast testing only
		if !*flagLongTest {
			btpParams.SlotsToCoeffsParameters.LogSlots = btpParams.BootstrappingParameters.LogN() - 1
			btpParams.CoeffsToSlotsParameters.LogSlots = btpParams.BootstrappingParameters.LogN() - 1

			// Corrects the message ratio to take into account the smaller number of slots and keep the same precision
			btpParams.Mod1ParametersLiteral.LogMessageRatio += 16 - params.LogN()
		}
		t.Logf("ParamsN2: LogN=%d/LogSlots=%d/LogQP=%f", params.LogN(), params.LogMaxSlots(), params.LogQP())

		// sk := rlwe.NewKeyGenerator(btpParams.BootstrappingParameters).GenSecretKeyNew()

		t.Log("Generating Bootstrapping Keys")
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

		require.NoError(t, err)

		prng, _ := sampling.NewKeyedPRNG([]byte{'a', 'b', 'c'})

		samplerQ := ring.NewUniformSampler(prng, ringQP.RingQ)
		samplerP := ring.NewUniformSampler(prng, ringQP.RingP)

		polyP := samplerP.ReadNew()
		polyQ := samplerQ.ReadNew()
		levelP := polyP.Level()
		fmt.Println("polyQ: ", polyQ.Level(), "polyP: ", polyP.Level())

		polyQ2 := ringQP.RingQ.NewPoly()
		num_run := 100
		warmup := 100

		// polyQ := plaintext.
		for i := 0; i < warmup; i++ {
			be.ModDownQPtoQ(polyQ.Level(), polyP.Level(), polyQ, polyP, polyQ2)
		}

		now := time.Now()
		for i := 0; i < num_run; i++ {
			be.ModDownQPtoQ(polyQ.Level(), polyP.Level(), polyQ, polyP, polyQ2)
		}
		fmt.Println("moddown: ", time.Since(now))

		moduli := append(ringQP.RingQ.ModuliChain(), ringQP.RingP.ModuliChain()...)
		ringQP_singlering, _ := ring.NewRing(ringQP.N(), moduli)
		samplerQP := ring.NewUniformSampler(prng, ringQP_singlering)

		polyQP2 := ringQP_singlering.NewPoly()
		polyQP3 := ringQP_singlering.NewPoly()
		rQPLvl := ringQP_singlering.Level()

		polyQP := samplerQP.ReadNew()
		var dur time.Duration

		// for i := 0; i < 10; i++ {
		// 	polyQP = samplerQP.ReadNew()
		// 	for idx, _ := range ringQP.RingP.ModuliChain() {
		// 		ringQP_singlering.AtLevel(rQPLvl-idx).DivRoundByLastModulus(polyQP, polyQP)
		// 		polyQP.Resize(rQPLvl - idx - 1)
		// 	}
		// }

		ring_now := ringQP_singlering.AtLevel(rQPLvl)
		ring_now.DivRoundByLastModulusAux(polyQP, polyQP3, levelP)
		polyQP3.Resize(rQPLvl - 1)
		ring_now.DivRoundByLastModulus(polyQP, polyQP2)
		polyQP2.Resize(rQPLvl - 1)

		// polyQP2.Resize(rQPLvl - 0 - 1)
		// polyQP3.Resize(rQPLvl - 0 - 1)
		// polyQP = samplerQP.ReadNew()

		// 	ring_now.DivRoundByLastModulus(polyQP, polyQP2)
		// 	ring_now.DivRoundByLastModulusAux(polyQP, polyQP3)
		diffs := make([]uint64, len(polyQP2.Coeffs[0]))
		for col_idx, v := range polyQP2.Coeffs[0] {
			diffs[col_idx] = v - polyQP3.Coeffs[0][col_idx]
		}

		for row_idx, row := range polyQP2.Coeffs {
			for col_idx, v := range row {
				if diffs[col_idx] != v-polyQP3.Coeffs[row_idx][col_idx] {
					fmt.Println("not equal at ", row_idx, col_idx)
					fmt.Println("v = ", v, "v2 = ", polyQP3.Coeffs[row_idx][col_idx], "diff = ", diffs[col_idx])
					break
				}
			}
		}
		fmt.Println(polyQP2.Coeffs[0][:5])
		fmt.Println(polyQP3.Coeffs[0][:5])

		// for row_idx, row := range polyQP3.Coeffs {

		// 	for col_idx, val := range row {
		// 		val2 := polyQP2.Coeffs[row_idx][col_idx]
		// 		if val != val2 {
		// 			fmt.Println("not equal")
		// 			break
		// 		}
		// 	}
		// }
		for idx, _ := range ringQP.RingP.ModuliChain() {
			ring_now := ringQP_singlering.AtLevel(rQPLvl - idx)
			for i := 0; i < warmup; i++ {

				ring_now.DivRoundByLastModulusAux(polyQP, polyQP2, levelP)
				// polyQP2.Resize(rQPLvl - idx - 1)

			}

			now = time.Now()
			for i := 0; i < num_run; i++ {
				ring_now.DivRoundByLastModulusAux(polyQP, polyQP2, levelP)
			}
			dur += time.Since(now)
			fmt.Println("rescale: ", idx, ",", time.Since(now))

			ringQP_singlering.AtLevel(rQPLvl-idx).DivRoundByLastModulus(polyQP, polyQP)
			polyQP.Resize(rQPLvl - idx - 1)
		}
		fmt.Println("rescale: ", dur)

	})
}

func TestMyRescaling(t *testing.T) {
	t.Run("RescalingLevelByLevel", func(t *testing.T) {
		schemeParamsLit := test_rescale_set1.SchemeParams
		btpParamsLit := ParametersLiteral{}

		if *flagLongTest {
			schemeParamsLit.LogN = 16
		}

		params, err := ckks.NewParametersFromLiteral(schemeParamsLit)
		require.Nil(t, err)

		btpParamsLit.LogN = utils.Pointy(params.LogN())

		btpParams, err := NewParametersFromLiteral(params, btpParamsLit)
		require.Nil(t, err)

		// Insecure params for fast testing only
		if !*flagLongTest {
			btpParams.SlotsToCoeffsParameters.LogSlots = btpParams.BootstrappingParameters.LogN() - 1
			btpParams.CoeffsToSlotsParameters.LogSlots = btpParams.BootstrappingParameters.LogN() - 1

			// Corrects the message ratio to take into account the smaller number of slots and keep the same precision
			btpParams.Mod1ParametersLiteral.LogMessageRatio += 16 - params.LogN()
		}
		t.Logf("ParamsN2: LogN=%d/LogSlots=%d/LogQP=%f", params.LogN(), params.LogMaxSlots(), params.LogQP())

		// sk := rlwe.NewKeyGenerator(btpParams.BootstrappingParameters).GenSecretKeyNew()

		// t.Log("Generating Bootstrapping Keys")
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

		require.NoError(t, err)

		prng, _ := sampling.NewKeyedPRNG([]byte{'a', 'b', 'c'})

		samplerQ := ring.NewUniformSampler(prng, ringQP.RingQ)
		samplerP := ring.NewUniformSampler(prng, ringQP.RingP)

		polyP := samplerP.ReadNew()
		polyQ := samplerQ.ReadNew()
		fmt.Println("polyQ: ", polyQ.Level(), "polyP: ", polyP.Level())

		polyQ2 := ringQP.RingQ.NewPoly()

		// polyQ := plaintext.

		warmup := 20
		test_round := 220
		fmt.Println("warmup=", warmup)
		fmt.Println("test_round=", test_round-warmup)

		var acc_time_modd time.Duration
		for i := 0; i < test_round; i++ {
			now := time.Now()
			be.ModDownQPtoQ(polyQ.Level(), polyP.Level(), polyQ, polyP, polyQ2)
			if i >= warmup {
				acc_time_modd += time.Since(now)
			}
		}
		fmt.Println("moddown: total time=", acc_time_modd, "average time=", acc_time_modd/(time.Duration(test_round-warmup)))

		moduli := append(ringQP.RingQ.ModuliChain(), ringQP.RingP.ModuliChain()...)
		ringQP_singlering, _ := ring.NewRing(ringQP.N(), moduli)
		samplerQP := ring.NewUniformSampler(prng, ringQP_singlering)

		// polyQP2 := ringQP_singlering.NewPoly()
		rQPLvl := ringQP_singlering.Level()

		var total_time_rescal time.Duration

		for i := 0; i < test_round; i++ {
			polyQP := samplerQP.ReadNew()

			for idx, _ := range ringQP.RingP.ModuliChain() {
				now := time.Now()
				ringQP_singlering.AtLevel(rQPLvl-idx).DivRoundByLastModulusAux(polyQP, polyQP, 0)
				if i >= warmup {
					total_time_rescal += time.Since(now)
				}

				polyQP.Resize(rQPLvl - idx - 1)
			}
		}
		fmt.Println("rescale: tatal time=", total_time_rescal, "total=", total_time_rescal/(time.Duration(test_round-warmup)))
	})
}

func TestBootstrapping(t *testing.T) {

	t.Run("BootstrappingWithoutRingDegreeSwitch", func(t *testing.T) {

		schemeParamsLit := test_rescale_set1.SchemeParams     // testPrec45
		btpParamsLit := test_rescale_set1.BootstrappingParams //ParametersLiteral{}

		if *flagLongTest {
			schemeParamsLit.LogN = 16
		}

		params, err := ckks.NewParametersFromLiteral(schemeParamsLit)
		require.Nil(t, err)

		btpParamsLit.LogN = utils.Pointy(params.LogN())

		btpParams, err := NewParametersFromLiteral(params, btpParamsLit)
		require.Nil(t, err)

		// Insecure params for fast testing only
		if !*flagLongTest {
			btpParams.SlotsToCoeffsParameters.LogSlots = btpParams.BootstrappingParameters.LogN() - 1
			btpParams.CoeffsToSlotsParameters.LogSlots = btpParams.BootstrappingParameters.LogN() - 1

			// Corrects the message ratio to take into account the smaller number of slots and keep the same precision
			btpParams.Mod1ParametersLiteral.LogMessageRatio += 16 - params.LogN()
		}

		t.Logf("ParamsN2: LogN=%d/LogSlots=%d/LogQP=%f", params.LogN(), params.LogMaxSlots(), params.LogQP())

		sk := rlwe.NewKeyGenerator(btpParams.BootstrappingParameters).GenSecretKeyNew()

		t.Log("Generating Bootstrapping Keys")
		btpKeys, _, err := btpParams.GenEvaluationKeys(sk)
		require.NoError(t, err)

		evaluator, err := NewEvaluator(btpParams, btpKeys)
		require.NoError(t, err)

		ecd := ckks.NewEncoder(params)
		enc := rlwe.NewEncryptor(params, sk)
		dec := rlwe.NewDecryptor(params, sk)

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

		t.Run("Bootstrapping", func(t *testing.T) {

			plaintext := ckks.NewPlaintext(params, 0)
			ecd.Encode(values, plaintext)

			ctQ0, err := enc.EncryptNew(plaintext)
			require.NoError(t, err)

			// Checks that the input ciphertext is at the level 0
			require.True(t, ctQ0.Level() == 0)

			// Bootstrapps the ciphertext
			now := time.Now()
			ctQL, err := evaluator.Bootstrap(ctQ0)
			fmt.Println("BTS time:", time.Since(now))
			require.NoError(t, err)

			// Checks that the output ciphertext is at the max level of paramsN1
			require.True(t, ctQL.Level() == params.MaxLevel())
			require.True(t, ctQL.Scale.Equal(params.DefaultScale()))

			verifyTestVectorsBootstrapping(params, ecd, dec, values, ctQL, t)
		})
	})

	// t.Run("BootstrappingWithRingDegreeSwitch", func(t *testing.T) {

	// 	schemeParamsLit := testPrec45
	// 	btpParamsLit := ParametersLiteral{}

	// 	if *flagLongTest {
	// 		schemeParamsLit.LogN = 16
	// 	}

	// 	schemeParamsLit.LogNthRoot = schemeParamsLit.LogN + 1
	// 	schemeParamsLit.LogN--

	// 	params, err := ckks.NewParametersFromLiteral(schemeParamsLit)
	// 	require.Nil(t, err)

	// 	btpParamsLit.LogN = utils.Pointy(params.LogN() + 1)

	// 	btpParams, err := NewParametersFromLiteral(params, btpParamsLit)
	// 	require.Nil(t, err)

	// 	// Insecure params for fast testing only
	// 	if !*flagLongTest {
	// 		btpParams.SlotsToCoeffsParameters.LogSlots = btpParams.BootstrappingParameters.LogN() - 1
	// 		btpParams.CoeffsToSlotsParameters.LogSlots = btpParams.BootstrappingParameters.LogN() - 1

	// 		// Corrects the message ratio to take into account the smaller number of slots and keep the same precision
	// 		btpParams.Mod1ParametersLiteral.LogMessageRatio += 16 - params.LogN()
	// 	}

	// 	t.Logf("Params: LogN=%d/LogSlots=%d/LogQP=%f", btpParams.ResidualParameters.LogN(), btpParams.ResidualParameters.LogMaxSlots(), btpParams.ResidualParameters.LogQP())
	// 	t.Logf("BTPParams: LogN=%d/LogSlots=%d/LogQP=%f", btpParams.BootstrappingParameters.LogN(), btpParams.BootstrappingParameters.LogMaxSlots(), btpParams.BootstrappingParameters.LogQP())

	// 	sk := rlwe.NewKeyGenerator(params).GenSecretKeyNew()

	// 	t.Log("Generating Bootstrapping Keys")
	// 	btpKeys, _, err := btpParams.GenEvaluationKeys(sk)
	// 	require.Nil(t, err)

	// 	evaluator, err := NewEvaluator(btpParams, btpKeys)
	// 	require.Nil(t, err)

	// 	ecd := ckks.NewEncoder(params)
	// 	enc := rlwe.NewEncryptor(params, sk)
	// 	dec := rlwe.NewDecryptor(params, sk)

	// 	values := make([]complex128, params.MaxSlots())
	// 	for i := range values {
	// 		values[i] = sampling.RandComplex128(-1, 1)
	// 	}

	// 	values[0] = complex(0.9238795325112867, 0.3826834323650898)
	// 	values[1] = complex(0.9238795325112867, 0.3826834323650898)
	// 	if len(values) > 2 {
	// 		values[2] = complex(0.9238795325112867, 0.3826834323650898)
	// 		values[3] = complex(0.9238795325112867, 0.3826834323650898)
	// 	}

	// 	t.Run("N1ToN2->Bootstrapping->N2ToN1", func(t *testing.T) {

	// 		plaintext := ckks.NewPlaintext(params, 0)
	// 		ecd.Encode(values, plaintext)

	// 		ctQ0, err := enc.EncryptNew(plaintext)
	// 		require.NoError(t, err)

	// 		// Checks that the input ciphertext is at the level 0
	// 		require.True(t, ctQ0.Level() == 0)

	// 		// Bootstrapps the ciphertext
	// 		ctQL, err := evaluator.Bootstrap(ctQ0)

	// 		if err != nil {
	// 			t.Fatal(err)
	// 		}

	// 		// Checks that the output ciphertext is at the max level of params
	// 		require.True(t, ctQL.Level() == params.MaxLevel())
	// 		require.True(t, ctQL.Scale.Equal(params.DefaultScale()))

	// 		verifyTestVectorsBootstrapping(params, ecd, dec, values, ctQL, t)

	// 	})
	// })

	// t.Run("BootstrappingPackedWithoutRingDegreeSwitch", func(t *testing.T) {

	// 	schemeParamsLit := testPrec45
	// 	btpParamsLit := ParametersLiteral{}

	// 	if *flagLongTest {
	// 		schemeParamsLit.LogN = 16
	// 	}

	// 	btpParamsLit.LogN = utils.Pointy(schemeParamsLit.LogN)

	// 	params, err := ckks.NewParametersFromLiteral(schemeParamsLit)
	// 	require.Nil(t, err)

	// 	// Insecure params for fast testing only
	// 	if !*flagLongTest {
	// 		// Corrects the message ratio to take into account the smaller number of slots and keep the same precision
	// 		btpParamsLit.LogMessageRatio = utils.Pointy(DefaultLogMessageRatio + (16 - params.LogN()))
	// 	}

	// 	sk := rlwe.NewKeyGenerator(params).GenSecretKeyNew()

	// 	ecd := ckks.NewEncoder(params)
	// 	enc := rlwe.NewEncryptor(params, sk)
	// 	dec := rlwe.NewDecryptor(params, sk)

	// 	for _, sparsity := range []int{0, 1, 2} {
	// 		btpParamsLit.LogSlots = utils.Pointy(*btpParamsLit.LogN - 1 - sparsity)
	// 		btpParams, err := NewParametersFromLiteral(params, btpParamsLit)
	// 		require.Nil(t, err)

	// 		t.Logf("Params: LogN=%d/LogSlots=%d/LogQP=%f", btpParams.ResidualParameters.LogN(), btpParams.ResidualParameters.LogMaxSlots(), btpParams.ResidualParameters.LogQP())
	// 		t.Logf("BTPParams: LogN=%d/LogSlots=%d/LogQP=%f", btpParams.BootstrappingParameters.LogN(), btpParams.BootstrappingParameters.LogMaxSlots(), btpParams.BootstrappingParameters.LogQP())

	// 		t.Log("Generating Bootstrapping Keys for LogSlots")
	// 		btpKeys, _, err := btpParams.GenEvaluationKeys(sk)
	// 		require.Nil(t, err)

	// 		evaluator, err := NewEvaluator(btpParams, btpKeys)
	// 		require.Nil(t, err)

	// 		for _, slotOffset := range []int{0, 1, 2, 3} {
	// 			logMaxSlots := params.LogMaxSlots() - slotOffset
	// 			logMaxSlots = utils.Min(logMaxSlots, btpParams.LogMaxSlots())
	// 			values := make([]complex128, 1<<logMaxSlots)
	// 			for i := range values {
	// 				values[i] = sampling.RandComplex128(-1, 1)
	// 			}

	// 			values[0] = complex(0.9238795325112867, 0.3826834323650898)
	// 			values[1] = complex(0.9238795325112867, 0.3826834323650898)
	// 			if len(values) > 2 {
	// 				values[2] = complex(0.9238795325112867, 0.3826834323650898)
	// 				values[3] = complex(0.9238795325112867, 0.3826834323650898)
	// 			}

	// 			pt := ckks.NewPlaintext(params, 0)
	// 			pt.LogDimensions = ring.Dimensions{Rows: 0, Cols: logMaxSlots}

	// 			cts := make([]rlwe.Ciphertext, 11)
	// 			for i := range cts {

	// 				require.NoError(t, ecd.Encode(utils.RotateSlice(values, i), pt))

	// 				ct, err := enc.EncryptNew(pt)
	// 				require.NoError(t, err)

	// 				cts[i] = *ct
	// 			}

	// 			if cts, err = evaluator.BootstrapMany(cts); err != nil {
	// 				t.Fatal(err)
	// 			}

	// 			for i := range cts {
	// 				// Checks that the output ciphertext is at the max level of paramsN1
	// 				require.True(t, cts[i].Level() == params.MaxLevel())
	// 				require.True(t, cts[i].Scale.Equal(params.DefaultScale()))

	// 				verifyTestVectorsBootstrapping(params, ecd, dec, utils.RotateSlice(values, i), &cts[i], t)
	// 			}
	// 		}
	// 	}
	// })

	// t.Run("BootstrappingPackedWithRingDegreeSwitch", func(t *testing.T) {

	// 	schemeParamsLit := testPrec45
	// 	btpParamsLit := ParametersLiteral{}

	// 	if *flagLongTest {
	// 		schemeParamsLit.LogN = 16
	// 	}

	// 	btpParamsLit.LogN = utils.Pointy(schemeParamsLit.LogN)
	// 	schemeParamsLit.LogNthRoot = schemeParamsLit.LogN + 1
	// 	schemeParamsLit.LogN -= 3

	// 	params, err := ckks.NewParametersFromLiteral(schemeParamsLit)
	// 	require.Nil(t, err)

	// 	// Insecure params for fast testing only
	// 	if !*flagLongTest {
	// 		// Corrects the message ratio to take into account the smaller number of slots and keep the same precision
	// 		btpParamsLit.LogMessageRatio = utils.Pointy(DefaultLogMessageRatio + (16 - params.LogN()))
	// 	}

	// 	sk := rlwe.NewKeyGenerator(params).GenSecretKeyNew()

	// 	ecd := ckks.NewEncoder(params)
	// 	enc := rlwe.NewEncryptor(params, sk)
	// 	dec := rlwe.NewDecryptor(params, sk)

	// 	for _, sparsity := range []int{0, 1, 2} {
	// 		btpParamsLit.LogSlots = utils.Pointy(*btpParamsLit.LogN - 1 - sparsity)
	// 		btpParams, err := NewParametersFromLiteral(params, btpParamsLit)
	// 		require.Nil(t, err)

	// 		t.Logf("Params: LogN=%d/LogSlots=%d/LogQP=%f", btpParams.ResidualParameters.LogN(), btpParams.ResidualParameters.LogMaxSlots(), btpParams.ResidualParameters.LogQP())
	// 		t.Logf("BTPParams: LogN=%d/LogSlots=%d/LogQP=%f", btpParams.BootstrappingParameters.LogN(), btpParams.BootstrappingParameters.LogMaxSlots(), btpParams.BootstrappingParameters.LogQP())

	// 		t.Log("Generating Bootstrapping Keys for LogSlots")
	// 		btpKeys, _, err := btpParams.GenEvaluationKeys(sk)
	// 		require.Nil(t, err)

	// 		evaluator, err := NewEvaluator(btpParams, btpKeys)
	// 		require.Nil(t, err)

	// 		for _, slotOffset := range []int{0, 1, 2, 3} {
	// 			logMaxSlots := params.LogMaxSlots() - slotOffset
	// 			logMaxSlots = utils.Min(logMaxSlots, btpParams.LogMaxSlots())
	// 			values := make([]complex128, 1<<logMaxSlots)
	// 			for i := range values {
	// 				values[i] = sampling.RandComplex128(-1, 1)
	// 			}

	// 			values[0] = complex(0.9238795325112867, 0.3826834323650898)
	// 			values[1] = complex(0.9238795325112867, 0.3826834323650898)
	// 			if len(values) > 2 {
	// 				values[2] = complex(0.9238795325112867, 0.3826834323650898)
	// 				values[3] = complex(0.9238795325112867, 0.3826834323650898)
	// 			}

	// 			pt := ckks.NewPlaintext(params, 0)
	// 			pt.LogDimensions = ring.Dimensions{Rows: 0, Cols: logMaxSlots}

	// 			cts := make([]rlwe.Ciphertext, 11)
	// 			for i := range cts {

	// 				require.NoError(t, ecd.Encode(utils.RotateSlice(values, i), pt))

	// 				ct, err := enc.EncryptNew(pt)
	// 				require.NoError(t, err)

	// 				cts[i] = *ct
	// 			}

	// 			if cts, err = evaluator.BootstrapMany(cts); err != nil {
	// 				t.Fatal(err)
	// 			}

	// 			for i := range cts {
	// 				// Checks that the output ciphertext is at the max level of paramsN1
	// 				require.True(t, cts[i].Level() == params.MaxLevel())
	// 				require.True(t, cts[i].Scale.Equal(params.DefaultScale()))

	// 				verifyTestVectorsBootstrapping(params, ecd, dec, utils.RotateSlice(values, i), &cts[i], t)
	// 			}
	// 		}
	// 	}
	// })

	// t.Run("BootstrappingWithRingTypeSwitch", func(t *testing.T) {

	// 	schemeParamsLit := testPrec45
	// 	schemeParamsLit.RingType = ring.ConjugateInvariant
	// 	btpParamsLit := ParametersLiteral{}

	// 	if *flagLongTest {
	// 		schemeParamsLit.LogN = 16
	// 	}

	// 	btpParamsLit.LogN = utils.Pointy(schemeParamsLit.LogN)
	// 	schemeParamsLit.LogNthRoot = schemeParamsLit.LogN + 1
	// 	schemeParamsLit.LogN--

	// 	params, err := ckks.NewParametersFromLiteral(schemeParamsLit)
	// 	require.Nil(t, err)

	// 	btpParams, err := NewParametersFromLiteral(params, btpParamsLit)
	// 	require.Nil(t, err)

	// 	// Insecure params for fast testing only
	// 	if !*flagLongTest {
	// 		// Corrects the message ratio to take into account the smaller number of slots and keep the same precision
	// 		btpParams.Mod1ParametersLiteral.LogMessageRatio += 16 - params.LogN()
	// 	}

	// 	t.Logf("Params: LogN=%d/LogSlots=%d/LogQP=%f", btpParams.ResidualParameters.LogN(), btpParams.ResidualParameters.LogMaxSlots(), btpParams.ResidualParameters.LogQP())
	// 	t.Logf("BTPParams: LogN=%d/LogSlots=%d/LogQP=%f", btpParams.BootstrappingParameters.LogN(), btpParams.BootstrappingParameters.LogMaxSlots(), btpParams.BootstrappingParameters.LogQP())

	// 	sk := rlwe.NewKeyGenerator(params).GenSecretKeyNew()

	// 	t.Log("Generating Bootstrapping Keys")
	// 	btpKeys, _, err := btpParams.GenEvaluationKeys(sk)
	// 	require.Nil(t, err)

	// 	evaluator, err := NewEvaluator(btpParams, btpKeys)
	// 	require.Nil(t, err)

	// 	ecd := ckks.NewEncoder(params)
	// 	enc := rlwe.NewEncryptor(params, sk)
	// 	dec := rlwe.NewDecryptor(params, sk)

	// 	values := make([]float64, params.MaxSlots())
	// 	for i := range values {
	// 		values[i] = sampling.RandFloat64(-1, 1)
	// 	}

	// 	values[0] = 0.9238795325112867
	// 	values[1] = 0.9238795325112867
	// 	if len(values) > 2 {
	// 		values[2] = 0.9238795325112867
	// 		values[3] = 0.9238795325112867
	// 	}

	// 	t.Run("ConjugateInvariant->Standard->Bootstrapping->Standard->ConjugateInvariant", func(t *testing.T) {

	// 		plaintext := ckks.NewPlaintext(params, 0)
	// 		require.NoError(t, ecd.Encode(values, plaintext))

	// 		ctLeftQ0, err := enc.EncryptNew(plaintext)
	// 		require.NoError(t, err)
	// 		ctRightQ0, err := enc.EncryptNew(plaintext)
	// 		require.NoError(t, err)

	// 		// Checks that the input ciphertext is at the level 0
	// 		require.True(t, ctLeftQ0.Level() == 0)
	// 		require.True(t, ctRightQ0.Level() == 0)

	// 		// Bootstraps the ciphertext
	// 		ctLeftQL, ctRightQL, err := evaluator.EvaluateConjugateInvariant(ctLeftQ0, ctRightQ0)

	// 		require.NoError(t, err)

	// 		// Checks that the output ciphertext is at the max level of paramsN1
	// 		require.True(t, ctLeftQL.Level() == params.MaxLevel())
	// 		require.True(t, ctLeftQL.Scale.Equal(params.DefaultScale()))

	// 		verifyTestVectorsBootstrapping(params, ecd, dec, values, ctLeftQL, t)

	// 		require.True(t, ctRightQL.Level() == params.MaxLevel())
	// 		require.True(t, ctRightQL.Scale.Equal(params.DefaultScale()))
	// 		verifyTestVectorsBootstrapping(params, ecd, dec, values, ctRightQL, t)
	// 	})
	// })
}

func verifyTestVectorsBootstrapping(params ckks.Parameters, encoder *ckks.Encoder, decryptor *rlwe.Decryptor, valuesWant, element interface{}, t *testing.T) {
	precStats := ckks.GetPrecisionStats(params, encoder, decryptor, valuesWant, element, 0, false)
	if *printPrecisionStats {
		t.Log(precStats.String())
	}

	rf64 := precStats.AVGLog2Prec.Real
	if64 := precStats.AVGLog2Prec.Imag

	minPrec := math.Log2(params.DefaultScale().Float64()) - float64(params.LogN()+2)
	if minPrec < 0 {
		minPrec = 0
	}

	minPrec -= 10

	require.GreaterOrEqual(t, rf64, minPrec)
	require.GreaterOrEqual(t, if64, minPrec)
}
