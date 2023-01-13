/*
 * Class NativeUtils is published under the The MIT License:
 *
 * Copyright (c) 2012 Adam Heinrich <adam@adamh.cz>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
package be.ugent.jgaborator.util;

import java.io.*;
import java.nio.file.FileSystemNotFoundException;
import java.nio.file.FileSystems;
import java.nio.file.Files;
import java.nio.file.ProviderNotFoundException;
import java.nio.file.StandardCopyOption;
import java.util.Arrays;

/**
 * A simple library class which helps with loading dynamic libraries stored in the
 * JAR archive. These libraries usually contain implementation of some methods in
 * native code (using JNI - Java Native Interface).
 * 
 * @see <a href="http://adamheinrich.com/blog/2012/how-to-load-native-jni-library-from-jar">http://adamheinrich.com/blog/2012/how-to-load-native-jni-library-from-jar</a>
 * @see <a href="https://github.com/adamheinrich/native-utils">https://github.com/adamheinrich/native-utils</a>
 *
 */
public class ZigNativeUtils {
 
    /**
     * The minimum length a prefix for a file has to have according to {@link File#createTempFile(String, String)}}.
     */
    private static final int MIN_PREFIX_LENGTH = 3;
    public static final String NATIVE_FOLDER_PATH_PREFIX = "zignativeutils";

    /**
     * Temporary directory which will contain the DLLs.
     */
    private static File temporaryDir;

    /**
     * Private constructor - this class will never be instanced
     */
    private ZigNativeUtils() {
    }

    /**
     *
     * @return
     * @see https://github.com/trustin/os-maven-plugin/blob/master/src/main/java/kr/motd/maven/os/Detector.java#L163
     * @see https://github.com/trustin/os-maven-plugin/blob/master/src/main/java/kr/motd/maven/os/Detector.java#L192
     */
    private static String ziggifyOS(){
        String os = System.getProperty("os.name").toLowerCase();
        if(os.startsWith("mac") || os.startsWith("osx")){
            os = "macos";
        }else if (os.contains("indows")){
            os = "windows";
        }else if (os.contains("linux")){
            os = "linux";
        }else if (os.startsWith("freebsd")) {
            os = "freebsd";
        } else  if (os.startsWith("openbsd")) {
            os = "openbsd";
        } else if  (os.startsWith("netbsd")) {
            os = "netbsd";
        } else {
            System.err.println("OS " + os + " not recognized will. Will use prefix '" + os + "_'" );
        }
        return os;
    }

    private static String ziggifyArch() {
        String value = System.getProperty("os.arch").toLowerCase();
        if (value.matches("^(x8664|amd64|ia32e|em64t|x64)$")) {
            return "x86_64";
        }
        if (value.matches("^(x8632|x86|i[3-6]86|ia32|x32)$")) {
            return "i386";
        }
        if (value.matches("^(ia64w?|itanium64)$")) {
            return "itanium_64";
        }
        if ("ia64n".equals(value)) {
            return "itanium_32";
        }
        if (value.matches("^(sparc|sparc32)$")) {
            return "sparc_32";
        }
        if (value.matches("^(sparcv9|sparc64)$")) {
            return "sparc_64";
        }
        if (value.matches("^(arm|arm32)$")) {
            return "arm_32";
        }
        if ("aarch64".equals(value)) {
            return "aarch64";
        }
        if (value.matches("^(mips|mips32)$")) {
            return "mips";
        }
        if (value.matches("^(mipsel|mips32el)$")) {
            return "mipsel";
        }
        if ("mips64".equals(value)) {
            return "mips64";
        }
        if ("mips64el".equals(value)) {
            return "mips64el";
        }
        if (value.matches("^(ppc|ppc32)$")) {
            return "powerpc";
        }
        if (value.matches("^(ppcle|ppc32le)$")) {
            return "ppcle_32";
        }
        if ("ppc64".equals(value)) {
            return "ppc_64";
        }
        if ("ppc64le".equals(value)) {
            return "ppcle_64";
        }
        if ("s390".equals(value)) {
            return "s390_32";
        }
        if ("s390x".equals(value)) {
            return "s390_64";
        }
        if (value.matches("^(riscv|riscv32)$")) {
            return "riscv";
        }
        if ("riscv64".equals(value)) {
            return "riscv64";
        }
        if ("e2k".equals(value)) {
            return "e2k";
        }
        if ("loongarch64".equals(value)) {
            return "loongarch_64";
        }

        return value;
    }

    private static String makeZigPrefix(){
        String os = ziggifyOS();
        String arch = ziggifyArch();

        if(!Arrays.asList(zigOperatingSystems).contains(os)){
            System.err.println("OS " + os + " not in the list of OSes recognized by Zig:");
            for(String zigOS : Arrays.asList(zigOperatingSystems)){
                System.err.println("\t" + zigOS);
            };
        }

        if(!Arrays.asList(zigArchitectures).contains(arch)){
            System.err.println("Architecture " + arch + " not in the list of architectures recognized by Zig:");
            for(String zigArch : Arrays.asList(zigArchitectures)){
                System.err.println("\t" + zigArch);
            };
        }

        return os + "_" + arch;
    }


    /**
     * Load a library with OS and architecture detection. The idea is that libraries are found
     * @param path
     * @throws IOException
     */
    public static void loadLibraryFromJarWithOSDetection(String path) throws IOException{
        String[] parts = path.split("/");
        String original_filename = (parts.length > 1) ? parts[parts.length - 1] : null;
        String prefixed_filename = makeZigPrefix()  + "_" + original_filename;
        String prefixed_path = path.replace(original_filename,prefixed_filename);

        System.err.println("Try to load JNI library " + prefixed_path + " from JAR archive.");
        ZigNativeUtils.loadLibraryFromJar(prefixed_path);
        System.err.println("Native library library " + prefixed_path + " loaded!");
    }


    /**
     * Loads library from current JAR archive
     * 
     * The file from JAR is copied into system temporary directory and then loaded. The temporary file is deleted after
     * exiting.
     * Method uses String as filename because the pathname is "abstract", not system-dependent.
     * 
     * @param path The path of file inside JAR as absolute path (beginning with '/'), e.g. /package/File.ext
     * @throws IOException If temporary file creation or read/write operation fails
     * @throws IllegalArgumentException If source file (param path) does not exist
     * @throws IllegalArgumentException If the path is not absolute or if the filename is shorter than three characters
     * (restriction of {@link File#createTempFile(java.lang.String, java.lang.String)}).
     * @throws FileNotFoundException If the file could not be found inside the JAR.
     */
    public static void loadLibraryFromJar(String path) throws IOException {
 
        if (null == path || !path.startsWith("/")) {
            throw new IllegalArgumentException("The path has to be absolute (start with '/').");
        }
 
        // Obtain filename from path
        String[] parts = path.split("/");
        String filename = (parts.length > 1) ? parts[parts.length - 1] : null;
 
        // Check if the filename is okay
        if (filename == null || filename.length() < MIN_PREFIX_LENGTH) {
            throw new IllegalArgumentException("The filename has to be at least 3 characters long.");
        }
 
        // Prepare temporary file
        if (temporaryDir == null) {
            temporaryDir = createTempDirectory(NATIVE_FOLDER_PATH_PREFIX);
            temporaryDir.deleteOnExit();
        }

        File temp = new File(temporaryDir, filename);

        try (InputStream is = ZigNativeUtils.class.getResourceAsStream(path)) {
            Files.copy(is, temp.toPath(), StandardCopyOption.REPLACE_EXISTING);
        } catch (IOException e) {
            temp.delete();
            throw e;
        } catch (NullPointerException e) {
            temp.delete();
            throw new FileNotFoundException("File " + path + " was not found inside JAR.");
        }

        try {
            System.load(temp.getAbsolutePath());
        } finally {
            if (isPosixCompliant()) {
                // Assume POSIX compliant file system, can be deleted after loading
                temp.delete();
            } else {
                // Assume non-POSIX, and don't delete until last file descriptor closed
                temp.deleteOnExit();
            }
        }
    }

    private static boolean isPosixCompliant() {
        try {
            return FileSystems.getDefault()
                    .supportedFileAttributeViews()
                    .contains("posix");
        } catch (FileSystemNotFoundException
                | ProviderNotFoundException
                | SecurityException e) {
            return false;
        }
    }

    private static File createTempDirectory(String prefix) throws IOException {
        String tempDir = System.getProperty("java.io.tmpdir");
        File generatedDir = new File(tempDir, prefix + System.nanoTime());
        
        if (!generatedDir.mkdir())
            throw new IOException("Failed to create temp directory " + generatedDir.getName());
        
        return generatedDir;
    }

    /**
     * A list of operating systems known by Zig.
     * Generate with command 'zig targets'
     * @see https://ziglang.org/
     */
    private final static String[] zigOperatingSystems = {"freestanding",
            "ananas",
            "cloudabi",
            "dragonfly",
            "freebsd",
            "fuchsia",
            "ios",
            "kfreebsd",
            "linux",
            "lv2",
            "macos",
            "netbsd",
            "openbsd",
            "solaris",
            "windows",
            "zos",
            "haiku",
            "minix",
            "rtems",
            "nacl",
            "aix",
            "cuda",
            "nvcl",
            "amdhsa",
            "ps4",
            "elfiamcu",
            "tvos",
            "watchos",
            "mesa3d",
            "contiki",
            "amdpal",
            "hermit",
            "hurd",
            "wasi",
            "emscripten",
            "uefi",
            "opencl",
            "glsl450",
            "vulkan",
            "plan9",
            "other"};
    /**
     * A list of architectures known by Zig.
     * Generate with command 'zig targets'
     * @see https://ziglang.org/
     */
    private static final String[] zigArchitectures = {"arm",
            "armeb",
            "aarch64",
            "aarch64_be",
            "aarch64_32",
            "arc",
            "avr",
            "bpfel",
            "bpfeb",
            "csky",
            "hexagon",
            "m68k",
            "mips",
            "mipsel",
            "mips64",
            "mips64el",
            "msp430",
            "powerpc",
            "powerpcle",
            "powerpc64",
            "powerpc64le",
            "r600",
            "amdgcn",
            "riscv32",
            "riscv64",
            "sparc",
            "sparcv9",
            "sparcel",
            "s390x",
            "tce",
            "tcele",
            "thumb",
            "thumbeb",
            "i386",
            "x86_64",
            "xcore",
            "nvptx",
            "nvptx64",
            "le32",
            "le64",
            "amdil",
            "amdil64",
            "hsail",
            "hsail64",
            "spir",
            "spir64",
            "kalimba",
            "shave",
            "lanai",
            "wasm32",
            "wasm64",
            "renderscript32",
            "renderscript64",
            "ve",
            "spu_2",
            "spirv32",
            "spirv64"
    };

}