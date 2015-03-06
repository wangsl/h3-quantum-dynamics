#!/bin/bash

# $Id$

function svnId()
{
    for f in $*; do
	svn propset svn:keywords "Id" $f
    done
}

svn add *.[FChm] *.pl *.sh *.inc *.make 

svnId *.[FChm] *.pl *.sh *.inc *.make 







