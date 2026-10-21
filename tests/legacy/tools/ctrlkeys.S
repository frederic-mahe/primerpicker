# CTRLKEYS.COM: a DOS TSR for the probes of PRIMER.exe. DOSBox-X's
# autotype cannot hold Ctrl down, so this program turns function keys
# into Ctrl+key codes in the BIOS keyboard buffer: F5 Ctrl+Home,
# F6 Ctrl+End, F7 Ctrl+PgUp, F8 Ctrl+PgDn, F9 Ctrl+Left, F10 Ctrl+Right.
#
# Build (GNU binutils): as --32 -o ctrlkeys.o ctrlkeys.S
#   ld -m elf_i386 -Ttext 0x100 --oformat binary -o CTRLKEYS.COM ctrlkeys.o
# (probe_picking.sh builds it). After the BIOS handler (INT 9) has
# stored a key, the handler below replaces the last key of the buffer
# when it is in the table (PRIMER.exe's own use of F5 to F10, if any,
# is then lost).
        .code16
        .text
        .globl _start
_start:
        jmp     install

old9:   .word   0, 0
table:  .word   0x3f00, 0x7700   # F5  -> Ctrl+Home
        .word   0x4000, 0x7500   # F6  -> Ctrl+End
        .word   0x4100, 0x8400   # F7  -> Ctrl+PgUp
        .word   0x4200, 0x7600   # F8  -> Ctrl+PgDn
        .word   0x4300, 0x7300   # F9  -> Ctrl+Left
        .word   0x4400, 0x7400   # F10 -> Ctrl+Right
table_end:

handler:
        pushf
        lcall   *%cs:old9
        push    %ax
        push    %bx
        push    %si
        push    %ds
        mov     $0x40, %ax
        mov     %ax, %ds
        mov     0x1c, %bx          # tail
        cmp     0x1a, %bx          # head
        je      done
        sub     $2, %bx
        cmp     0x80, %bx          # buffer start
        jae     have
        mov     0x82, %bx          # buffer end
        sub     $2, %bx
have:   mov     (%bx), %ax
        mov     $table, %si
next:   cmp     $table_end, %si
        jae     done
        cmp     %cs:(%si), %ax
        je      found
        add     $4, %si
        jmp     next
found:  mov     %cs:2(%si), %ax
        mov     %ax, (%bx)
done:   pop     %ds
        pop     %si
        pop     %bx
        pop     %ax
        iret

install:
        mov     $0x3509, %ax
        int     $0x21
        mov     %bx, old9
        mov     %es, old9+2
        mov     $0x2509, %ax
        mov     $handler, %dx
        int     $0x21
        mov     $install, %dx
        add     $15, %dx
        shr     $4, %dx
        mov     $0x3100, %ax
        int     $0x21
