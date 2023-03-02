load("@rules_cc//cc:defs.bzl", "cc_binary")

cc_binary(
    name = "lm",
    srcs = [
        "main.cpp",
        ],
    includes = [
      "src",
    ],
    deps = [
        "//src/mylib:mylib",
    ],
    visibility = ["//visibility:public"],
)
