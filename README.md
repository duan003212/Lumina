
# Lumina: High-Performance Monte Carlo Ray Tracer

![C++ Standard](https://img.shields.io/badge/C%2B%2B-17-blue.svg?style=flat&logo=c%2B%2B)
![License](https://img.shields.io/badge/License-MIT-green.svg)
![Build Status](https://img.shields.io/badge/build-passing-brightgreen.svg)

**Lumina** 是一个基于 **C++17** 开发的高性能物理渲染引擎（Physically Based Renderer）。它不依赖任何图形 API（如 OpenGL/DirectX），而是基于**蒙特卡洛路径追踪算法（Monte Carlo Path Tracing）**，从零开始实现了光线与物体交互的数学模型。

本项目旨在探索计算机图形学的底层原理，并利用现代 C++ 特性与并行计算技术（OpenMP）优化渲染效率。

## 🖼️ Gallery (渲染展示)

![1114673c-af33-4004-ae7a-d1792c8e0a6d](https://github.com/user-attachments/assets/035c48e8-e0fa-4621-81e6-ee5f30b537c3)


## ✨ Key Features (核心特性)

*   **物理真实感渲染 (PBR):**
    *   基于 **Global Illumination** (全局光照) 模拟光线的多次弹射。
    *   支持 **Lambertian** (漫反射)、**Metal** (金属反射，支持模糊度)、**Dielectric** (电介质/玻璃，基于 Snell's Law 和 Schlick 近似) 等多种材质。
*   **高性能架构:**
    *   **Multi-threading:** 利用 **OpenMP** 实现多线程并行渲染，在多核 CPU 上实现近乎线性的加速比。
    *   **Modern C++:** 广泛使用 `std::shared_ptr` 进行内存管理，利用 `std::optional` 和 `lambda` 表达式增强代码鲁棒性。
*   **高级相机模型:**
    *   支持 **Anti-Aliasing (MSAA)** 多重采样抗锯齿。
    *   模拟薄透镜相机，支持 **Depth of Field (景深)** 和 **Defocus Blur (散焦模糊)** 效果。
*   **数学与算法:**
    *   从零构建向量数学库 (`vec3`)。
    *   实现 **Russian Roulette** (俄罗斯轮盘赌) 策略优化路径终止条件。

## 🛠️ Tech Stack (技术栈)

*   **Language:** C++17
*   **Parallel Computing:** OpenMP
*   **Build System:** CMake (Version 3.10+) / Make
*   **Output Format:** PPM (Portable Pixel Map)

## 🚀 Getting Started (如何运行)

### Prerequisites (环境要求)
*   支持 C++17 的编译器 (GCC 7+, Clang 5+, MSVC 19.14+)
*   CMake (可选)
  
### Build with G++ (快速测试)

直接通过命令行编译（确保开启 O3 优化和 OpenMP）：

```bash
g++ -std=c++17 -O3 -fopenmp main.cpp -o raytracer
./raytracer > image.ppm
```

## 📊 Performance Benchmarks (性能基准)

测试环境: Intel Core i7-13620H (10 Cores), 16GB RAM, Image Size: 1200x800, Samples: 100.

| Mode | Render Time | Speedup |
| :--- | :--- | :--- |
| Single Thread | 86.4s | 1.0x |
| **OpenMP (8 Threads)** | **11.2s** | **7.7x** |

*注：并行优化显著消除了渲染瓶颈，利用率达到 CPU 峰值的 95% 以上。*


## 📚 References (参考资料)

本项目深受 Peter Shirley 的经典著作启发：
*   [_Ray Tracing in One Weekend_](https://raytracing.github.io/books/RayTracingInOneWeekend.html)

---
