#!/bin/bash
aclocal -I m4 --install
autoconf --force
automake --force --copy --add-missing
libtoolize --copy --install --force
