#!/bin/bash

# $Id$

function svnId()
{
    for f in $*; do
	svn propset svn:keywords "Id" $f
    done
}

svn add *.[FChm] *.pl *.sh Makefile-BKMP2Mex.make

svnId add *.[FChm] *.pl *.sh Makefile-BKMP2Mex.make

