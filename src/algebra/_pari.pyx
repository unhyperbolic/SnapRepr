#
#
# This file initializes pari on load (_pari_init) and provides _pari_eval.
# This file is supposed to be accesed through pari.py
#
# _pari_eval is a wrapper to gp_read_str
# The input and output is a python string (using GENtostr).
#
# The output is what the gp command line would give if the input is typed in.

# For details, look at 
# \url{http://pari.math.u-bordeaux.fr/pub/pari/manuals/2.3.3/libpari.pdf} 

cdef extern from "stdlib.h":
    void free(void*)

cdef extern from "setjmp.h":
    struct __jump_buf_tag:
        pass
    ctypedef __jump_buf_tag jmp_buf
    int setjmp(jmp_buf)

#cdef extern from "cheated_pari_headers/pari.h":
cdef extern from "pari/pari.h":
    ctypedef unsigned long pari_sp
    ctypedef long *GEN
    ctypedef unsigned long size_t
    void pari_init(size_t, unsigned long)
    void pari_init_opts(size_t, unsigned long, unsigned long)
    GEN gp_read_str(char *)
    char* GENtostr(GEN)
    cdef extern pari_sp avma

#cdef extern from "cheated_pari_headers/paripriv.h":
cdef extern from "pari/paripriv.h":
    cdef enum:
        INIT_DFTm=4
        INIT_SIGm=2
    ctypedef struct gp_data:
        jmp_buf env
        pass
    cdef extern gp_data *GP_DATA



def _pari_eval(str py_input):
    cdef GEN pari_result
    cdef char* cstr_result

    cdef pari_sp av

    global avma
    global GP_DATA

    av = avma
#    if setjmp(GP_DATA.env):
#        avma = av
#        raise Exception, "PARI ERROR in %s" % py_input

    bytes_input = bytes(py_input)
    
    pari_result = gp_read_str(bytes_input)
    cstr_result = GENtostr(pari_result)

    avma = av

    cdef bytes py_string = cstr_result
    free(cstr_result)
    return py_string

def _pari_init():
    pari_init_opts(5000000,0000000,INIT_DFTm)

_pari_init()
