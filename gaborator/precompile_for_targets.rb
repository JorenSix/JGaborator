platforms = [
"x86_64-linux-gnu",
"x86_64-windows-gnu",
"x86_64-macos-gnu",
"aarch64-macos-gnu",
"aarch64-linux-gnu",
"aarch64-linux-musl",
"aarch64-windows-gnu",
"i386-linux-gnu",
"i386-windows-gnu",
"riscv64-linux-musl",
]

extensions = {"windows" => "dll","macos" => "dylib","linux"=> "so"}
basename = "libjgaborator"
output_folder = "../src/main/resources/jni/"

platforms.each_with_index do |platform , index|
    os = platform.split("-")[1]
	extension = extensions[os]
	ouput_file = File.join(output_folder,"#{"%02d" % index}_#{platform}_#{basename}.#{extension}")
	puts ouput_file
	system("export ZIG_OUTPUT_FILE=#{ouput_file} && export ZIG_TARGET_PLATFORM=#{platform} && make zig")
end