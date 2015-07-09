prefix=$(realpath $(dir $(lastword $(MAKEFILE_LIST))))
export prefix
export atimeincludedir=$(prefix)/atime/include

build:
	$(MAKE) -C dgcpm build
	$(MAKE) -C atime build
	$(MAKE) -C kp build

install:
	- mkdir $(prefix)/lib $(prefix)/include $(prefix)/bin
	$(MAKE) -C dgcpm install
	$(MAKE) -C atime install
	$(MAKE) -C kp install

clean:
	$(MAKE) -C dgcpm clean
	$(MAKE) -C atime clean
	$(MAKE) -C kp clean

uninstall:
	$(MAKE) -C dgcpm uninstall
	$(MAKE) -C atime uninstall
	$(MAKE) -C kp uninstall
	- rmdir $(prefix)/lib $(prefix)/include $(prefix)/bin
