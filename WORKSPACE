workspace(name = "lm")

load("@bazel_tools//tools/build_defs/repo:http.bzl", "http_archive")

# Setup rules_foreign_cc (CMake integration)
http_archive(
    name = "rules_foreign_cc",
    strip_prefix = "rules_foreign_cc-8d540605805fb69e24c6bf5dc885b0403d74746a", # 0.9.0
    url = "https://github.com/bazelbuild/rules_foreign_cc/archive/8d540605805fb69e24c6bf5dc885b0403d74746a.tar.gz",
)

load("@rules_foreign_cc//foreign_cc:repositories.bzl", "rules_foreign_cc_dependencies")

rules_foreign_cc_dependencies()

http_archive(
      name = "name",
      build_file = "//:third_party/name.BUILD",
      sha256 = "8df0079cdbe179748a18d44731af62a245a45ebf5085223dc03133954c662973",
      strip_prefix = "name-prefix",
      urls = ["https://github.com/....tar.gz"],
)

