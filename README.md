# Prerequisites
GCC is required for measuring CPU cycles. The following setting is for go to enable CGO:
```sh
$ go env -w CGO_ENABLED 1
```

# Correctness
<!-- The benchmarking for section 5.1 are performed in `circuits/ckks/bootstrapping` using the following codes:
```sh
$ go test -v -benchmark=BenchmarkBootstrap -run=^$
``` -->

To perform the correctness test, run the following commands in the same directory:
```sh
$ go test -v -run TestCorrectness
```
It tests the results between `ModDown` and consecutive `Rescale`. There are small noises introduced between these two operations.


# CPU cycles
The codes for measuring cpu cycles are in `examples/cpucycles`. To run it, type `go run main.go` in the shell.


# Test in Boostrapping
To test the running time of Moddown algorithms in Section 5.2, change to directory `circuits/ckks/bootstrapping` and run
```sh
 go test -run TestMyRescaling
```
The parameters are set to `L=24,k=5` by default, one can change the value of `L` by modifying the `LogQ` (line 132 in `circuits/ckks/bootstrapping/default_parameters.go`).