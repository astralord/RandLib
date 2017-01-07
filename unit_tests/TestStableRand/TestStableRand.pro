#-------------------------------------------------
#
# Project created by QtCreator 2017-01-06T22:36:45
#
#-------------------------------------------------

QT       += testlib

QT       -= gui

TARGET = tst_testpdf
CONFIG   += console
CONFIG   -= app_bundle

TEMPLATE = app


SOURCES += \
    tst_stable.cpp
DEFINES += SRCDIR=\\\"$$PWD/\\\"

win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../../bin/RandLib/release/ -lrandlib
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../../bin/RandLib/debug/ -lrandlib
else:unix: LIBS += -L$$PWD/../../../bin/ -lrandlib

INCLUDEPATH += $$PWD/../../../RandLib
DEPENDPATH += $$PWD/../../../RandLib

CONFIG += c++11

HEADERS += \
    data.h
