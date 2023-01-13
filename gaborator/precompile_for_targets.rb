platforms = ["aarch64-linux-musl","aarch64-windows-gnu","aarch64-macos-gnu","i386-linux-gnu","i386-windows-gnu","riscv64-linux-musl","x86_64-linux-musl","x86_64-windows-gnu","x86_64-macos-gnu"]

extensions = {"windows" => "dll","macos" => "dylib","linux"=> "so"}
basename = "libjgaborator"
output_folder = "../src/main/resources/jni/"

platforms.each do |platform|
	
	os = platform.split("-")[1]
	arch = platform.split("-")[0]
	extension = extensions[os]

	ouput_file = File.join(output_folder,"#{os}_#{arch}_#{basename}.#{extension}")
	puts ouput_file
	system("export ZIG_OUTPUT_FILE=#{ouput_file} && export ZIG_TARGET_PLATFORM=#{platform} && make zig")
end