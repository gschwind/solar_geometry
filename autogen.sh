#!/bin/bash
mkdir -p m4
aclocal -I m4 --install
libtoolize --copy --install --force
autoconf --force
automake --force --copy --add-missing
