package bootstrapping

import (
	"fmt"
	"testing"
	"math"

	"github.com/stretchr/testify/assert"

	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	"github.com/tuneinsight/lattigo/v6/utils"
	"github.com/tuneinsight/lattigo/v6/utils/sampling"
)

func TestCorrectness(t *testing.T) {
	schemeParamsLit := N16QP1767H32768H32.SchemeParams
	btpParamsLit := ParametersLiteral{}

	if *flagLongTest {
		schemeParamsLit.LogN = 16
	}

	params, err := ckks.NewParametersFromLiteral(schemeParamsLit)
	if err != nil {
		panic(err)
	}

	btpParamsLit.LogN = utils.Pointy(params.LogN())

	btpParams, err := NewParametersFromLiteral(params, btpParamsLit)

	// Insecure params for fast testing only
	if !*flagLongTest {
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

	// ctQ0, err := enc.EncryptNew(plaintext)

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
	// polyQP2 := ringQP_singlering.NewPoly()
	// polyQP3 := ringQP_singlering.NewPoly()
	polyQP := samplerQP.ReadNew()
	rQPLvl := ringQP_singlering.Level()
	// num_run := 100
	// warmup := 100

	for pi, poly_mod_pi := range polyQ.Coeffs {
		for j, v := range poly_mod_pi {
			polyQP.Coeffs[pi][j] = v
		}
	}
	len_p := len(polyQ.Coeffs)
	for qi, poly_mod_qi := range polyP.Coeffs {
		for j, v := range poly_mod_qi {
			polyQP.Coeffs[len_p+qi][j] = v
		}
	}

	be.ModDownQPtoQ(polyQ.Level(), polyP.Level(), polyQ, polyP, polyQ2)

	for idx, _ := range ringQP.RingP.ModuliChain() {
		ringQP_singlering.AtLevel(rQPLvl-idx).DivRoundByLastModulusAux(polyQP, polyQP, 0)
	}

	// there are noises between the results due to the final implementations
	// the noises are small
	for pi, poly_mod_pi := range polyQ2.Coeffs {
		for j, v := range poly_mod_pi {
			assert.Less(t, float64(v)-float64(polyQP.Coeffs[pi][j]), 1024.0)
			assert.Less(t, float64(polyQP.Coeffs[pi][j])-float64(v), 1024.0)
		}
	}


	// furthermore, the difference is in R_Q. In RNS, we just check that every column has a same difference
	diff := make([]uint64, polyQ.N())
	for j, v := range polyQ2.Coeffs[0] {
		diff[j] = polyQP.Coeffs[0][j] - v
	}

	for i, poly_mod_pi := range polyQ2.Coeffs[1:] {
		for j, v := range poly_mod_pi {
			assert.Equal(t, diff[j], polyQP.Coeffs[i+1][j]-v)
		}
	}
	
	fmt.Println("former values:")
	fmt.Println(polyQ2.Coeffs[0][:10])
	fmt.Println(polyQP.Coeffs[0][:10])

	max_noise := 0.0
	for j, v := range polyQ2.Coeffs[0] {
		max_noise = max(max_noise, math.Abs(float64(polyQP.Coeffs[0][j])-float64(v)))
	}
	fmt.Println("max noise: ", max_noise)
}

