// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#pragma once

#define V_MAYORVERSION		1
#define V_MINORVERSION		0
#define V_BUGFIXVERSION		0
//#define V_BUILDVERSION		0

#define V_COMPANYNAME		"Lobo Lab\0"

#define V_FILEDESCRIPTION	"Viewer\0"

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

#define V_PRODUCTVERSION	"v" TOSTRING(V_MAYORVERSION) "." TOSTRING(V_MINORVERSION) "." TOSTRING(V_BUGFIXVERSION)
#define V_LEGALCOPYRIGHT	"Lobo Lab (lobolab.umbc.edu). All rights reserved.\0"
#define V_ORIGINALFILENAME	"viewer.exe\0"
#define V_PRODUCTNAME		"Evolution Viewer\0"
