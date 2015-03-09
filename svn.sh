#!/bin/bash

# $Id$

function svnId()
{
    for f in $*; do
	svn propset svn:keywords "Id" $f
    done
}

svn add *.[FChm] *.pl *.sh *.inc Makefile

svnId *.[FChm] *.pl *.sh *.inc Makefile

exit








