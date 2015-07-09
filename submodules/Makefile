prefix=$(realpath $(dir $(lastword $(MAKEFILE_LIST))))
export prefix

build:
	$(MAKE) -C dgcpm build

install:
	- mkdir $(prefix)/lib $(prefix)/include $(prefix)/bin
	$(MAKE) -C dgcpm install

clean:
	$(MAKE) -C dgcpm clean

uninstall:
	$(MAKE) -C dgcpm uninstall
	- rmdir $(prefix)/lib $(prefix)/include $(prefix)/bin
