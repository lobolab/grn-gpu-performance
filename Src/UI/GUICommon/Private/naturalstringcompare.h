/*
 * This software was written by people from OnShore Consulting services LLC
 * <info@sabgroup.com> and placed in the public domain.
 *
 * We reserve no legal rights to any of this. You are free to do
 * whatever you want with it. And we make no guarantee or accept
 * any claims on damages as a result of this.
 *
 * If you change the software, please help us and others improve the
 * code by sending your modifications to us. If you choose to do so,
 * your changes will be included under this license, and we will add
 * your name to the list of contributors.
*/

#pragma once

#include <QString>
#include <QStringList>

namespace NaturalStringCompare {

int naturalStringCompare( const QString & lhs, const QString & rhs, Qt::CaseSensitivity caseSensitive=Qt::CaseSensitive );
QStringList naturalStringSort( const QStringList & list, Qt::CaseSensitivity caseSensitive=Qt::CaseSensitive );
bool naturalStringCompareLessThan( const QString & lhs, const QString & rhs );
bool naturalStringCaseInsensitiveCompareLessThan( const QString & lhs, const QString & rhs );

} // namespace NaturalStringCompare
