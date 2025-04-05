# CPPFLAGS += -I htslib/
CFLAGS   += -g -Wall -O2  -std=c99
LDFLAGS  += $(LIBS) -lz -lm -lpthread
BUILD_DIR = build

BINARY = cornetto
OBJ = $(BUILD_DIR)/main.o \
      $(BUILD_DIR)/cornetto.o \
      $(BUILD_DIR)/depth_main.o \
      $(BUILD_DIR)/fixdir_main.o \
      $(BUILD_DIR)/boringbits_main.o \
	  $(BUILD_DIR)/bigenough_main.o \
      $(BUILD_DIR)/thread.o \
	  $(BUILD_DIR)/misc.o \
	  $(BUILD_DIR)/misc_p.o \
	  $(BUILD_DIR)/error.o \
      $(BUILD_DIR)/find_telomere.o \
      $(BUILD_DIR)/telomere_windows.o \
      $(BUILD_DIR)/telomere_breaks.o \
	  $(BUILD_DIR)/dotter.o \
	  $(BUILD_DIR)/paf.o \
	  $(BUILD_DIR)/sdict.o \
	  $(BUILD_DIR)/sdust.o \
	  $(BUILD_DIR)/assbed.o \
	  $(BUILD_DIR)/seq.o

ifdef asan
	CFLAGS += -fsanitize=address -fno-omit-frame-pointer
	LDFLAGS += -fsanitize=address -fno-omit-frame-pointer
endif

.PHONY: clean distclean test install uninstall

$(BINARY): $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) $(LDFLAGS) -o $@

$(BUILD_DIR)/main.o: src/main.c
	$(CC) $(CFLAGS) $(CPPFLAGS) $< -c -o $@

$(BUILD_DIR)/cornetto.o: src/cornetto.c src/misc.h src/error.h src/cornetto.h
	$(CC) $(CFLAGS) $(CPPFLAGS) $< -c -o $@

$(BUILD_DIR)/depth_main.o: src/depth_main.c src/error.h
	$(CC) $(CFLAGS) $(CPPFLAGS) $< -c -o $@

$(BUILD_DIR)/fixdir_main.o: src/fixdir_main.c src/khash.h src/error.h
	$(CC) $(CFLAGS) $(CPPFLAGS) $< -c -o $@

$(BUILD_DIR)/boringbits_main.o: src/boringbits_main.c src/error.h
	$(CC) $(CFLAGS) $(CPPFLAGS) $< -c -o $@

$(BUILD_DIR)/bigenough_main.o: src/bigenough_main.c src/error.h
	$(CC) $(CFLAGS) $(CPPFLAGS) $< -c -o $@

$(BUILD_DIR)/assbed.o: src/assbed.c src/error.h src/kseq.h
	$(CC) $(CFLAGS) $(CPPFLAGS) $< -c -o $@

$(BUILD_DIR)/seq.o: src/seq.c src/kseq.h
	$(CC) $(CFLAGS) $(CPPFLAGS) $< -c -o $@

$(BUILD_DIR)/thread.o: src/thread.c src/cornetto.h
	$(CC) $(CFLAGS) $(CPPFLAGS) $< -c -o $@

$(BUILD_DIR)/misc.o: src/misc.c src/misc.h
	$(CC) $(CFLAGS) $(CPPFLAGS) $< -c -o $@

$(BUILD_DIR)/misc_p.o: src/misc_p.c src/misc.h
	$(CC) $(CFLAGS) $(CPPFLAGS) $< -c -o $@

$(BUILD_DIR)/error.o: src/error.c src/error.h
	$(CC) $(CFLAGS) $(CPPFLAGS) $< -c -o $@

# temolere stuff
$(BUILD_DIR)/find_telomere.o: src/find_telomere.c
	$(CC) $(CFLAGS) $(CPPFLAGS) $< -c -o $@

$(BUILD_DIR)/telomere_windows.o: src/telomere_windows.c
	$(CC) $(CFLAGS) $(CPPFLAGS) $< -c -o $@

$(BUILD_DIR)/telomere_breaks.o: src/telomere_breaks.c
	$(CC) $(CFLAGS) $(CPPFLAGS) $< -c -o $@

# minidot
$(BUILD_DIR)/dotter.o: src/minidot/dotter.c src/minidot/eps.h  src/minidot/kvec.h  src/minidot/paf.h  src/minidot/sdict.h
	$(CC) $(CFLAGS) $(CPPFLAGS) $< -c -o $@

$(BUILD_DIR)/paf.o: src/minidot/paf.c src/minidot/kseq.h  src/minidot/paf.h
	$(CC) $(CFLAGS) $(CPPFLAGS) $< -c -o $@

$(BUILD_DIR)/sdict.o: src/minidot/sdict.c src/minidot/sdict.h  src/minidot/khash.h
	$(CC) $(CFLAGS) $(CPPFLAGS) $< -c -o $@

# sdust
$(BUILD_DIR)/sdust.o: src/sdust/sdust.c src/sdust/sdust.h src/sdust/kvec.h src/sdust/kdq.h src/sdust/kseq.h src/sdust/kalloc.h src/sdust/kseq.h src/sdust/ketopt.h
	$(CC) $(CFLAGS) $(CPPFLAGS) $< -c -o $@

# htslib/libhts.a:
# 	@if test -e $(BUILD_DIR)/lib/libhts.a; then \
# 		echo "htslib found at htslib/libhts.a"; \
# 	else \
# 		echo "htslib not found at htslib/libhts.a"; \
# 		echo "Please run 'scripts/install-hts.sh' first"; \
# 		exit 1; \
# 	fi

clean:
	rm -rf $(BINARY) $(BUILD_DIR)/*.o

# Delete all gitignored files (but not directories)
distclean: clean
	git clean -f -X
	rm -rf $(BUILD_DIR)/*

test: $(BINARY)
	./test/test.sh

