// Copyright (c) 2010-2020, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

#ifndef MFEM_SYCL_HPP
#define MFEM_SYCL_HPP

#include "../config/config.hpp"

#ifdef MFEM_USE_SYCL

#include "../general/debug.hpp"
#define MFEM_DEVICE
#define MFEM_LAMBDA
#define MFEM_HOST_DEVICE
#define MFEM_DEVICE_SYNC
#define MFEM_STREAM_SYNC
#define MFEM_GPU_CHECK(x)

// Define the SYCL inner threading macros
#define SYCL_SHARED
#define SYCL_SYNC_THREAD itm.barrier(sycl::access::fence_space::local_space);
#define SYCL_THREAD_ID(k) itm.get_local_id(k);
#define SYCL_THREAD_SIZE(k) itm.get_local_range(k);
#define SYCL_FOREACH_THREAD(i,k,N) \
    for(int i=itm.get_local_id(k); i<N; i+=itm.get_local_range(k))

#include "mem_manager.hpp"
#include "device.hpp"
#include "../linalg/dtensor.hpp"

namespace mfem
{

/// Return the default sycl::queue used by MFEM.
sycl::queue &SyclQueue();

#define SYCL_KERNEL(...) { \
    sycl::queue &Q = SyclQueue();\
    Q.submit([&](sycl::handler &h) {__VA_ARGS__}); \
    Q.wait();\
}

#define SYCL_FORALL(i,N,...) \
    ForallWrap1D(N, h, [=] (int i) {__VA_ARGS__})

/// The forall kernel body wrapper
template <typename BODY>
inline void ForallWrap1D(const int N, sycl::handler &h, BODY &&body)
{
   h.parallel_for(sycl::range<1>(N), [=](sycl::id<1> k) {body(k);});
   return;
}

////////////////////////////////////////////////////////////////////////////
// SYCL_FORALL with a 3D CUDA block, sycl::h_item<3> &itm
#define SYCL_FORALL_3D(i,N,X,Y,Z,...) \
   ForallWrap3D<X,Y,Z>(N, h, [=] (const int i, sycl::nd_item<3> itm) {__VA_ARGS__})

/// The forall kernel body wrapper
template <int X, int Y, int Z, typename BODY>
inline void ForallWrap3D(const int N, sycl::handler &h, BODY &&body)
{
   if (N == 0) { return; }
   constexpr int B = X*Y*Z;
   const int L = static_cast<int>(ceil(cbrt(N)));
   dbg("N:%d, B:%d, L:%d", N, B, L);
   const sycl::range<3> GRID(L*X,L*Y,L*Z);
   const sycl::range<3> BLCK(X,Y,Z);

   h.parallel_for(sycl::nd_range<3>(GRID, BLCK), [=](sycl::nd_item<3> itm)
   {
      SyKernel3D(N, body, itm);
   });
   return;
}

template <typename BODY> static
void SyKernel3D(const int N, BODY body, sycl::nd_item<3> itm)
{
   const int k = itm.get_group_linear_id();
   if (k >= N) { return; }
   body(k, itm);
}

/// Get the number of SYCL devices
int SyGetDeviceCount();

/** @brief Function that determines if an SYCL kernel should be used, based on
    the current mfem::Device configuration. */
inline bool DeviceCanUseSycl() { return Device::Allows(Backend::SYCL); }

} // namespace mfem

#endif // MFEM_USE_SYCL

#endif // MFEM_SYCL_HPP
