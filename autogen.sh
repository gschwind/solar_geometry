#!/bin/bash
aclocal -I m4 --install
libtoolize --copy --install --force
autoconf --force
automake --force --copy --add-missing
