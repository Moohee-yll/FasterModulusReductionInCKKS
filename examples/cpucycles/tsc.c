#include "tsc.h"
static uint64_t tsc() {
	uint64_t x;
	asm volatile (
		"mfence\n"  // Memory barrier
		"lfence\n"  // Load barrier
		"rdtsc\n"   // Read time-stamp counter
		"shl $32, %%rdx\n"  // Shift EDX left by 32 bits
		"or %%rdx, %%rax"   // Combine EDX and EAX into RAX
		: "=a" (x)  // Output: EAX (lower part of RAX) mapped to x
		:           // No inputs
		: "rdx"     // Clobbers EDX
	);
	return x;
}