#***********************************************************************
#*                   GNU Lesser General Public License
#*
#* This file is part of the GFDL Flexible Modeling System (FMS).
#*
#* FMS is free software: you can redistribute it and/or modify it under
#* the terms of the GNU Lesser General Public License as published by
#* the Free Software Foundation, either version 3 of the License, or (at
#* your option) any later version.
#*
#* FMS is distributed in the hope that it will be useful, but WITHOUT
#* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#* FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#* for more details.
#*
#* You should have received a copy of the GNU Lesser General Public
#* License along with FMS.  If not, see <http://www.gnu.org/licenses/>.
#***********************************************************************

# This is an automake file for the platform directory of the FMS
# package.

# Ed Hartnett 2/22/19

# Include .h files.
AM_CPPFLAGS = -I${top_builddir}/.mod -I${top_srcdir}/include

# Build this uninstalled convenience library.
noinst_LTLIBRARIES = libplatform.la

# The convenience library depends on its source.
libplatform_la_SOURCES = platform.F90

# Mod file depends on its o file, is built and then installed.
platform_mod.mod: platform.lo

BUILT_SOURCES = platform_mod.mod
include_HEADERS = platform_mod.mod

SUFFIXES = .mod _mod.mod
.F90.mod:
	$(PPFCCOMPILE) -c $<
	@cp ${builddir}/$@ ${top_builddir}/.mod/
.F90_mod.mod:
	$(PPFCCOMPILE) -c $<
	@cp ${builddir}/$@ ${top_builddir}/.mod/

CLEANFILES = *.mod
