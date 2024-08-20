#! /bin/bash
# ----------------------------------------------------------------------------

aclocal \
&& autoreconf --force --install --verbose 
autoreconf -i
