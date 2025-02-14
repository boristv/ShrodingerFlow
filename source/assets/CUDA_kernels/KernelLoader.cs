using System;
using System.IO;
using ManagedCuda;
using ManagedCuda.NVRTC;

namespace source.assets.CUDA_kernels
{
    public class KernelLoader
    {
        private static CudaContext ctx = new CudaContext(0);

        public static CudaKernel load_kernel(String kernelName)
        {
            byte[] ptx = prepare_kernel(kernelName);
            
            return ctx.LoadKernelPTX(ptx, kernelName);
        }
        
        private static byte[] prepare_kernel(string kernelName)
        {
            string fileToCompile = File.ReadAllText("./assets/CUDA_kernels/" + kernelName + ".cu");
            
            CudaRuntimeCompiler rtc = new CudaRuntimeCompiler(fileToCompile, kernelName);

            rtc.Compile(new string[]{});
            
            string log = rtc.GetLogAsString();

            Console.WriteLine(log);

            byte[] ptx = rtc.GetPTX();

            rtc.Dispose();

            return ptx;
        }
    }
}