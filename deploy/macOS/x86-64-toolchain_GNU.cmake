set(CMAKE_SYSTEM_NAME Generic)
set(CMAKE_SYSTEM_PROCESSOR X86)

set(TOOLCHAIN_NAME x86_64-apple-darwin22)
set(COMPILER_PATH "/gcc/bin")
set(TOOLCHAIN_PATH "/cctools/bin")
set(LLVM_PATH "/usr/lib/llvm-14/bin")

set(CMAKE_C_COMPILER "${COMPILER_PATH}/${TOOLCHAIN_NAME}-gcc")
set(CMAKE_CXX_COMPILER "${COMPILER_PATH}/${TOOLCHAIN_NAME}-g++")
set(CMAKE_ASM_COMPILER "${COMPILER_PATH}/${TOOLCHAIN_NAME}-gcc")
set(CMAKE_RC_COMPILER "${LLVM_PATH}/llvm-rc")
set(CMAKE_LINKER "${TOOLCHAIN_PATH}/${TOOLCHAIN_NAME}-ld")
set(CMAKE_ADDR2LINE "${LLVM_PATH}/llvm-addr2line")
set(CMAKE_AR "${COMPILER_PATH}/${TOOLCHAIN_NAME}-gcc-ar")
set(CMAKE_DLLTOOL "${LLVM_PATH}/llvm-dlltool")
set(CMAKE_MT "${LLVM_PATH}/llvm-mt")
set(CMAKE_NM "${COMPILER_PATH}/${TOOLCHAIN_NAME}-gcc-nm")
set(CMAKE_OBJCOPY "${LLVM_PATH}/llvm-objcopy")
set(CMAKE_OBJDUMP "${LLVM_PATH}/llvm-objdump")
set(CMAKE_RANLIB "${COMPILER_PATH}/${TOOLCHAIN_NAME}-gcc-ranlib")
set(CMAKE_READELF "${LLVM_PATH}/llvm-readelf")
set(CMAKE_STRIP "${TOOLCHAIN_PATH}/${TOOLCHAIN_NAME}-strip")
