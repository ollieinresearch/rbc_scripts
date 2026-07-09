import dedalus.core.evaluator as evaluator

def patched_get_file(self, **kw):
    exists = self.current_file.exists() if self.comm.rank == 0 else None
    exists = self.comm.bcast(exists, root=0)
    if not exists:
        self.create_current_file()
    return self.open_file(**kw)

evaluator.H5FileHandlerBase.get_file = patched_get_file

import runpy, sys
sys.argv = sys.argv[1:]
runpy.run_path(sys.argv[0], run_name="__main__")

