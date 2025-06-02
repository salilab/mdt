"""configure test to determine size of void *"""


def _check(context):
    context.Message("Determining sizeof(void *)... ")
    text = """
#include <stdio.h>

int main(void) {
  printf("%d", (int)sizeof(void *));
  return 0;
}
"""
    ret = context.TryRun(text, ".c")
    if ret[0] == 0:
        context.env.Exit("Could not run sizeof check program")

    # Workaround for dumb systems (e.g. wine) which insert stuff into stdout:
    result = ret[1].split()[-1]

    context.Result(result)
    return int(result)


def configure_check(env):
    custom_tests = {'CheckSizeof': _check}
    conf = env.Configure(custom_tests=custom_tests)
    env['SIZEOF_POINTER'] = conf.CheckSizeof()
    conf.Finish()
