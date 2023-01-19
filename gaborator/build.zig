const std = @import("std");

pub fn build(b: *std.build.Builder) void {
    const cflags = [_][]const u8{ "-Wall", "-Wextra", "-Werror=return-type", "-std=gnu11", "-O3", "-ffast-math" };
    const cppflags = [_][]const u8{ "-I/Users/joren/Downloads/jdk-17.0.5+8/include", "-Igaborator-1.7", "-Ipffft", "-I/Users/joren/Downloads/jdk-17.0.5+8/include/win32", "-std=c++11", "-DGABORATOR_USE_PFFFT", "-O3", "-ffast-math" };
    const lib = b.addSharedLibrary("jgaborator", null, b.version(0, 7, 0));
    lib.addCSourceFile("pffft/pffft.c", &cflags);
    lib.addCSourceFile("pffft/fftpack.c", &cflags);

    lib.addCSourceFile("jgaborator.cc", &cppflags);
    lib.linkLibC();
    lib.linkLibCpp();

    const target = b.standardTargetOptions(.{});
    const mode = b.standardReleaseOptions();

    lib.setTarget(target);

    lib.setBuildMode(mode);
    lib.install();
}
