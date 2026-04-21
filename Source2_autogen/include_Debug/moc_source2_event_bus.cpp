/****************************************************************************
** Meta object code from reading C++ file 'source2_event_bus.h'
**
** Created by: The Qt Meta Object Compiler version 69 (Qt 6.10.0)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../../source2_event_bus.h"
#include <QtCore/qmetatype.h>

#include <QtCore/qtmochelpers.h>

#include <memory>


#include <QtCore/qxptype_traits.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'source2_event_bus.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 69
#error "This file was generated using the moc from 6.10.0. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

#ifndef Q_CONSTINIT
#define Q_CONSTINIT
#endif

QT_WARNING_PUSH
QT_WARNING_DISABLE_DEPRECATED
QT_WARNING_DISABLE_GCC("-Wuseless-cast")
namespace {
struct qt_meta_tag_ZN12UQFFEventBusE_t {};
} // unnamed namespace

template <> constexpr inline auto UQFFEventBus::qt_create_metaobjectdata<qt_meta_tag_ZN12UQFFEventBusE_t>()
{
    namespace QMC = QtMocConstants;
    QtMocHelpers::StringRefStorage qt_stringData {
        "UQFFEventBus",
        "logMessage",
        "",
        "LogEntry",
        "entry",
        "systemParametersChanged",
        "systemName",
        "QJsonObject",
        "params",
        "computationCompleted",
        "source",
        "equationName",
        "results",
        "simulationFrameUpdated",
        "frame",
        "time",
        "fieldState",
        "validationResult",
        "testName",
        "passed",
        "error",
        "details",
        "apiDataReceived",
        "query",
        "apiSource",
        "data",
        "csvGenerated",
        "filepath",
        "rowCount",
        "broadcastMessage",
        "channel",
        "message",
        "computationRequested",
        "target",
        "method"
    };

    QtMocHelpers::UintData qt_methods {
        // Signal 'logMessage'
        QtMocHelpers::SignalData<void(const LogEntry &)>(1, 2, QMC::AccessPublic, QMetaType::Void, {{
            { 0x80000000 | 3, 4 },
        }}),
        // Signal 'systemParametersChanged'
        QtMocHelpers::SignalData<void(const QString &, const QJsonObject &)>(5, 2, QMC::AccessPublic, QMetaType::Void, {{
            { QMetaType::QString, 6 }, { 0x80000000 | 7, 8 },
        }}),
        // Signal 'computationCompleted'
        QtMocHelpers::SignalData<void(const QString &, const QString &, const QJsonObject &)>(9, 2, QMC::AccessPublic, QMetaType::Void, {{
            { QMetaType::QString, 10 }, { QMetaType::QString, 11 }, { 0x80000000 | 7, 12 },
        }}),
        // Signal 'simulationFrameUpdated'
        QtMocHelpers::SignalData<void(int, double, const QJsonObject &)>(13, 2, QMC::AccessPublic, QMetaType::Void, {{
            { QMetaType::Int, 14 }, { QMetaType::Double, 15 }, { 0x80000000 | 7, 16 },
        }}),
        // Signal 'validationResult'
        QtMocHelpers::SignalData<void(const QString &, bool, double, const QString &)>(17, 2, QMC::AccessPublic, QMetaType::Void, {{
            { QMetaType::QString, 18 }, { QMetaType::Bool, 19 }, { QMetaType::Double, 20 }, { QMetaType::QString, 21 },
        }}),
        // Signal 'apiDataReceived'
        QtMocHelpers::SignalData<void(const QString &, const QString &, const QJsonObject &)>(22, 2, QMC::AccessPublic, QMetaType::Void, {{
            { QMetaType::QString, 23 }, { QMetaType::QString, 24 }, { 0x80000000 | 7, 25 },
        }}),
        // Signal 'csvGenerated'
        QtMocHelpers::SignalData<void(const QString &, int)>(26, 2, QMC::AccessPublic, QMetaType::Void, {{
            { QMetaType::QString, 27 }, { QMetaType::Int, 28 },
        }}),
        // Signal 'broadcastMessage'
        QtMocHelpers::SignalData<void(const QString &, const QJsonObject &)>(29, 2, QMC::AccessPublic, QMetaType::Void, {{
            { QMetaType::QString, 30 }, { 0x80000000 | 7, 31 },
        }}),
        // Signal 'computationRequested'
        QtMocHelpers::SignalData<void(const QString &, const QString &, const QJsonObject &)>(32, 2, QMC::AccessPublic, QMetaType::Void, {{
            { QMetaType::QString, 33 }, { QMetaType::QString, 34 }, { 0x80000000 | 7, 8 },
        }}),
    };
    QtMocHelpers::UintData qt_properties {
    };
    QtMocHelpers::UintData qt_enums {
    };
    return QtMocHelpers::metaObjectData<UQFFEventBus, qt_meta_tag_ZN12UQFFEventBusE_t>(QMC::MetaObjectFlag{}, qt_stringData,
            qt_methods, qt_properties, qt_enums);
}
Q_CONSTINIT const QMetaObject UQFFEventBus::staticMetaObject = { {
    QMetaObject::SuperData::link<QObject::staticMetaObject>(),
    qt_staticMetaObjectStaticContent<qt_meta_tag_ZN12UQFFEventBusE_t>.stringdata,
    qt_staticMetaObjectStaticContent<qt_meta_tag_ZN12UQFFEventBusE_t>.data,
    qt_static_metacall,
    nullptr,
    qt_staticMetaObjectRelocatingContent<qt_meta_tag_ZN12UQFFEventBusE_t>.metaTypes,
    nullptr
} };

void UQFFEventBus::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    auto *_t = static_cast<UQFFEventBus *>(_o);
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: _t->logMessage((*reinterpret_cast<std::add_pointer_t<LogEntry>>(_a[1]))); break;
        case 1: _t->systemParametersChanged((*reinterpret_cast<std::add_pointer_t<QString>>(_a[1])),(*reinterpret_cast<std::add_pointer_t<QJsonObject>>(_a[2]))); break;
        case 2: _t->computationCompleted((*reinterpret_cast<std::add_pointer_t<QString>>(_a[1])),(*reinterpret_cast<std::add_pointer_t<QString>>(_a[2])),(*reinterpret_cast<std::add_pointer_t<QJsonObject>>(_a[3]))); break;
        case 3: _t->simulationFrameUpdated((*reinterpret_cast<std::add_pointer_t<int>>(_a[1])),(*reinterpret_cast<std::add_pointer_t<double>>(_a[2])),(*reinterpret_cast<std::add_pointer_t<QJsonObject>>(_a[3]))); break;
        case 4: _t->validationResult((*reinterpret_cast<std::add_pointer_t<QString>>(_a[1])),(*reinterpret_cast<std::add_pointer_t<bool>>(_a[2])),(*reinterpret_cast<std::add_pointer_t<double>>(_a[3])),(*reinterpret_cast<std::add_pointer_t<QString>>(_a[4]))); break;
        case 5: _t->apiDataReceived((*reinterpret_cast<std::add_pointer_t<QString>>(_a[1])),(*reinterpret_cast<std::add_pointer_t<QString>>(_a[2])),(*reinterpret_cast<std::add_pointer_t<QJsonObject>>(_a[3]))); break;
        case 6: _t->csvGenerated((*reinterpret_cast<std::add_pointer_t<QString>>(_a[1])),(*reinterpret_cast<std::add_pointer_t<int>>(_a[2]))); break;
        case 7: _t->broadcastMessage((*reinterpret_cast<std::add_pointer_t<QString>>(_a[1])),(*reinterpret_cast<std::add_pointer_t<QJsonObject>>(_a[2]))); break;
        case 8: _t->computationRequested((*reinterpret_cast<std::add_pointer_t<QString>>(_a[1])),(*reinterpret_cast<std::add_pointer_t<QString>>(_a[2])),(*reinterpret_cast<std::add_pointer_t<QJsonObject>>(_a[3]))); break;
        default: ;
        }
    }
    if (_c == QMetaObject::IndexOfMethod) {
        if (QtMocHelpers::indexOfMethod<void (UQFFEventBus::*)(const LogEntry & )>(_a, &UQFFEventBus::logMessage, 0))
            return;
        if (QtMocHelpers::indexOfMethod<void (UQFFEventBus::*)(const QString & , const QJsonObject & )>(_a, &UQFFEventBus::systemParametersChanged, 1))
            return;
        if (QtMocHelpers::indexOfMethod<void (UQFFEventBus::*)(const QString & , const QString & , const QJsonObject & )>(_a, &UQFFEventBus::computationCompleted, 2))
            return;
        if (QtMocHelpers::indexOfMethod<void (UQFFEventBus::*)(int , double , const QJsonObject & )>(_a, &UQFFEventBus::simulationFrameUpdated, 3))
            return;
        if (QtMocHelpers::indexOfMethod<void (UQFFEventBus::*)(const QString & , bool , double , const QString & )>(_a, &UQFFEventBus::validationResult, 4))
            return;
        if (QtMocHelpers::indexOfMethod<void (UQFFEventBus::*)(const QString & , const QString & , const QJsonObject & )>(_a, &UQFFEventBus::apiDataReceived, 5))
            return;
        if (QtMocHelpers::indexOfMethod<void (UQFFEventBus::*)(const QString & , int )>(_a, &UQFFEventBus::csvGenerated, 6))
            return;
        if (QtMocHelpers::indexOfMethod<void (UQFFEventBus::*)(const QString & , const QJsonObject & )>(_a, &UQFFEventBus::broadcastMessage, 7))
            return;
        if (QtMocHelpers::indexOfMethod<void (UQFFEventBus::*)(const QString & , const QString & , const QJsonObject & )>(_a, &UQFFEventBus::computationRequested, 8))
            return;
    }
}

const QMetaObject *UQFFEventBus::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *UQFFEventBus::qt_metacast(const char *_clname)
{
    if (!_clname) return nullptr;
    if (!strcmp(_clname, qt_staticMetaObjectStaticContent<qt_meta_tag_ZN12UQFFEventBusE_t>.strings))
        return static_cast<void*>(this);
    return QObject::qt_metacast(_clname);
}

int UQFFEventBus::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QObject::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 9)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 9;
    }
    if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 9)
            *reinterpret_cast<QMetaType *>(_a[0]) = QMetaType();
        _id -= 9;
    }
    return _id;
}

// SIGNAL 0
void UQFFEventBus::logMessage(const LogEntry & _t1)
{
    QMetaObject::activate<void>(this, &staticMetaObject, 0, nullptr, _t1);
}

// SIGNAL 1
void UQFFEventBus::systemParametersChanged(const QString & _t1, const QJsonObject & _t2)
{
    QMetaObject::activate<void>(this, &staticMetaObject, 1, nullptr, _t1, _t2);
}

// SIGNAL 2
void UQFFEventBus::computationCompleted(const QString & _t1, const QString & _t2, const QJsonObject & _t3)
{
    QMetaObject::activate<void>(this, &staticMetaObject, 2, nullptr, _t1, _t2, _t3);
}

// SIGNAL 3
void UQFFEventBus::simulationFrameUpdated(int _t1, double _t2, const QJsonObject & _t3)
{
    QMetaObject::activate<void>(this, &staticMetaObject, 3, nullptr, _t1, _t2, _t3);
}

// SIGNAL 4
void UQFFEventBus::validationResult(const QString & _t1, bool _t2, double _t3, const QString & _t4)
{
    QMetaObject::activate<void>(this, &staticMetaObject, 4, nullptr, _t1, _t2, _t3, _t4);
}

// SIGNAL 5
void UQFFEventBus::apiDataReceived(const QString & _t1, const QString & _t2, const QJsonObject & _t3)
{
    QMetaObject::activate<void>(this, &staticMetaObject, 5, nullptr, _t1, _t2, _t3);
}

// SIGNAL 6
void UQFFEventBus::csvGenerated(const QString & _t1, int _t2)
{
    QMetaObject::activate<void>(this, &staticMetaObject, 6, nullptr, _t1, _t2);
}

// SIGNAL 7
void UQFFEventBus::broadcastMessage(const QString & _t1, const QJsonObject & _t2)
{
    QMetaObject::activate<void>(this, &staticMetaObject, 7, nullptr, _t1, _t2);
}

// SIGNAL 8
void UQFFEventBus::computationRequested(const QString & _t1, const QString & _t2, const QJsonObject & _t3)
{
    QMetaObject::activate<void>(this, &staticMetaObject, 8, nullptr, _t1, _t2, _t3);
}
QT_WARNING_POP
