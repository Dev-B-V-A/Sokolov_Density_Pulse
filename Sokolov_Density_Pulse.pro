QT += core
QT -= gui

CONFIG += c++11

TARGET = Sokolov_Density_Pulse
CONFIG += console
CONFIG -= app_bundle

TEMPLATE = app

SOURCES += main.cpp \
    gas_params.cpp \
    solver.cpp \
    scheme.cpp \
    norma.cpp \
    report.cpp \
    init.cpp

HEADERS += \
    gas_params.h \
    solver.h \
    scheme.h \
    norma.h \
    report.h \
    init.h
