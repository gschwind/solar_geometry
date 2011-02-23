#!/bin/bash
aclocal
libtoolize --force
autoconf
automake --add-missing
