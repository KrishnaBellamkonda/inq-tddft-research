// Minimal per-device free-memory probe. Uses cudaMemGetInfo, which does NOT go
// through NVML — so it works even though nvidia-smi is broken on this box
// (driver/library version mismatch). Prints one line per device:
//   GPU <i> free_MB <free> total_MB <total>
#include <cstdio>
#include <cuda_runtime.h>

int main() {
  int n = 0;
  if (cudaGetDeviceCount(&n) != cudaSuccess) {
    std::printf("ERR cudaGetDeviceCount\n");
    return 1;
  }
  for (int i = 0; i < n; i++) {
    if (cudaSetDevice(i) != cudaSuccess) {
      std::printf("GPU %d ERR setDevice\n", i);
      continue;
    }
    size_t fr = 0, tot = 0;
    if (cudaMemGetInfo(&fr, &tot) != cudaSuccess) {
      std::printf("GPU %d ERR memInfo\n", i);
      continue;
    }
    std::printf("GPU %d free_MB %zu total_MB %zu\n", i,
                fr / (1024 * 1024), tot / (1024 * 1024));
  }
  return 0;
}
