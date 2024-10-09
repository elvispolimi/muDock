#pragma once

#include <cuda_runtime.h>
#include <mudock/log.hpp>
#include <stdexcept>

// TODO check -Wterminate
#define MUDOCK_CHECK(call) call

#define MUDOCK_CHECK_KERNELCALL() printf("no kernel check\n")
