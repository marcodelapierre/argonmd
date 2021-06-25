.PHONY: all
all: serial_one serial_multi


.PHONY: serial_one
serial_one:
	cd serial_one_source && $(MAKE) exe

.PHONY: serial_multi
serial_multi:
	cd serial_multi_sources && $(MAKE) exe


.PHONY: clean_obj_all
clean_obj_all:
	rm -f */*.o

.PHONY: clean_all
clean_all:
	rm -f */*.x */*.o
