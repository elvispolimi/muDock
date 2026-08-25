When I started porting the code to GPU, I noticed something strange. The naive Alpaka GPU version was actually faster than the highly optimized Native CUDA version. This was weird because the Alpaka version was still very basic. After looking into it, I discovered that the difference was just the random number generator (RNG). When I used the same Philox generator in both versions, the performance matched exactly as expected.

Later, when I tested the CPU versions, I found another problem. The Alpaka CPU version was about 28% slower than the normal Native C++ version. Since I disabled loop unrolling, the code structure was almost identical. So, I investigated and found two main reasons for this performance drop:

### 1. The Generator Type

Alpaka uses the `Philox` generator. This is great for GPUs, but the math is a bit slow on CPUs. The Native C++ CPU code uses `std::mt19937` (Mersenne Twister), which is natively very fast on CPU architectures.

### 2. The Missing Reference

Just changing the algorithm to Mersenne Twister in Alpaka did not fix the problem. I realized there was a big memory issue. In the Alpaka kernel, the code copied the RNG state by value from global memory:

```cpp
alpaka_rand_state l_state = state[global_thread_id];
```

For the Philox generator on the GPU, this state is only 16 bytes. Copying it is very fast. However, the Mersenne Twister state is bigger . Copying a lot of data for every single thread at every kernel launch was too heavy for the CPU memory and slowed everything down.

---

## The Solution

To fix both problems without breaking the single codebase rule (Write-Once), I used some C++17 features. I created a `rng_selector` that automatically chooses the right generator and the right way to pass the memory depending on the hardware.

```cpp
template<typename TAcc, typename Enable = void>
struct rng_selector {
    using type = alpaka::rand::Philox4x32x10;
    using local_ref = type; // Pass by value for GPU

};

template<typename TAcc>
struct rng_selector<TAcc, std::enable_if_t<std::is_same_v<alpaka::Platform<TAcc>, alpaka::PlatformCpu>>> {
    using type = alpaka::rand::engine::cpu::MersenneTwister;
    using local_ref = type&; // Pass by reference for CPU 

    };
```

I also added a simple `if constexpr` to make sure the code only saves the state back to memory if it is not a reference:

```cpp
alpaka_rand_local l_state = state[global_thread_id];
// ... math operations ...
if constexpr (!std::is_reference_v<alpaka_rand_local>) {
    state[global_thread_id] = l_state;
}
```

## Final Results

This fix worked perfectly. By using the right generator and stopping the memory copy, the CPU backends now access the state directly via reference.

The results were great. The 28% gap was completely closed. In fact, the Alpaka Serial backend is now **38% faster** than the original Native C++ Serial code (about 15,600 Evals/s compared to 11,300 Evals/s on the 100x1000 dataset). The better structure required for the GPU code (like block-stride loops) helped the CPU a lot too.
