#!/usr/bin/env python3
"""Clear the executable bit on an ELF's PT_GNU_STACK segment.

Newer WSL2 kernels refuse to load shared objects that request an executable
stack (Baysor's bundled libopenlibm.so triggers this). This flips that one flag
off in place - equivalent to `execstack -c` / `patchelf --clear-execstack`.
"""
import struct
import sys

PT_GNU_STACK = 0x6474E551

for path in sys.argv[1:]:
    with open(path, "r+b") as f:
        data = bytearray(f.read())
        assert data[:4] == b"\x7fELF", f"{path}: not ELF"
        is64 = data[4] == 2
        le = data[5] == 1
        end = "<" if le else ">"
        if is64:
            e_phoff = struct.unpack_from(end + "Q", data, 0x20)[0]
            e_phentsize = struct.unpack_from(end + "H", data, 0x36)[0]
            e_phnum = struct.unpack_from(end + "H", data, 0x38)[0]
            flags_off = 4  # p_flags immediately after p_type in Elf64_Phdr
        else:
            e_phoff = struct.unpack_from(end + "I", data, 0x1C)[0]
            e_phentsize = struct.unpack_from(end + "H", data, 0x2A)[0]
            e_phnum = struct.unpack_from(end + "H", data, 0x2C)[0]
            flags_off = 24  # p_flags is last field in Elf32_Phdr

        changed = False
        for i in range(e_phnum):
            ph = e_phoff + i * e_phentsize
            p_type = struct.unpack_from(end + "I", data, ph)[0]
            if p_type == PT_GNU_STACK:
                fo = ph + flags_off
                p_flags = struct.unpack_from(end + "I", data, fo)[0]
                if p_flags & 0x1:
                    struct.pack_into(end + "I", data, fo, p_flags & ~0x1)
                    changed = True
                    print(f"{path}: PT_GNU_STACK flags {p_flags:#x} -> {p_flags & ~0x1:#x}")
        if changed:
            f.seek(0)
            f.write(data)
        else:
            print(f"{path}: nothing to change")
