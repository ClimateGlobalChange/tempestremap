///////////////////////////////////////////////////////////////////////////////
///
///	\file    FunctionTimer.cpp
///	\author  Paul Ullrich
///	\version November 2, 2021
///
///	<remarks>
///		Copyright 2021 Paul Ullrich
///
///		This file is distributed as part of the Tempest source code package.
///		Permission is granted to use, copy, modify and distribute this
///		source code and its documentation under the terms of the GNU General
///		Public License.  This software is provided "as is" without express
///		or implied warranty.
///	</remarks>

#include "FunctionTimer.h"
#include "Exception.h"

#include <iostream>
#include <chrono>

///////////////////////////////////////////////////////////////////////////////

FunctionTimer::GroupDataMap FunctionTimer::m_mapGroupData;

///////////////////////////////////////////////////////////////////////////////

FunctionTimer::FunctionTimer(const char *szGroup) {

	// Start the timer
	m_fStopped = false;

	// Assign group name
	if (szGroup == NULL) {
		m_strGroup = "";
	} else {
		m_strGroup = szGroup;
	}

	// Assign start time
	m_tpStartTime = std::chrono::high_resolution_clock::now();
}

///////////////////////////////////////////////////////////////////////////////

void FunctionTimer::Reset() {
	m_fStopped = false;
	m_tpStartTime = std::chrono::high_resolution_clock::now();
}

///////////////////////////////////////////////////////////////////////////////

std::chrono::microseconds FunctionTimer::Time(bool fDone) {

	if (!m_fStopped) {
		m_tpStopTime = std::chrono::high_resolution_clock::now();
	}

	std::chrono::microseconds msTime =
		std::chrono::duration_cast<std::chrono::microseconds>(
			m_tpStopTime - m_tpStartTime);

	// If no name associated with this timer, ignore fDone.
	if (m_strGroup == "") {
		return msTime;
	}

	// Add the time to the group record
	if ((fDone) && (!m_fStopped)) {
		m_fStopped = true;

		GroupDataMap::iterator iter;

		iter = m_mapGroupData.find(m_strGroup);

		// Add to existing group record
		if (iter != m_mapGroupData.end()) {
			iter->second.iTotalTime += msTime;
			iter->second.nEntries++;

		// Create new group record
		} else {
			GroupDataPair gdp;
			gdp.first = m_strGroup;
			gdp.second.iTotalTime = msTime;
			gdp.second.nEntries = 1;

			m_mapGroupData.insert(gdp);
		}
	}

	return msTime;
}

///////////////////////////////////////////////////////////////////////////////

std::chrono::microseconds FunctionTimer::StopTime() {
	return Time(true);
}

///////////////////////////////////////////////////////////////////////////////

const FunctionTimer::TimerGroupData & FunctionTimer::GetGroupTimeRecord(
	const char *szName
) {
	GroupDataMap::iterator iter;

	iter = m_mapGroupData.find(szName);

	// Retrieve existing group record
	if (iter != m_mapGroupData.end()) {
		return iter->second;

	// Group record does not exist
	} else {
		_EXCEPTION1("Group time record %s does not exist.", szName);
	}
}

///////////////////////////////////////////////////////////////////////////////

std::chrono::microseconds FunctionTimer::GetTotalGroupTime(const char *szName) {

	GroupDataMap::iterator iter;

	iter = m_mapGroupData.find(szName);

	// Retrieve existing group record
	if (iter != m_mapGroupData.end()) {
		const TimerGroupData & tgd = iter->second;

		return (tgd.iTotalTime);

	// Group record does not exist
	} else {
		return std::chrono::microseconds(0);
	}
}

///////////////////////////////////////////////////////////////////////////////

std::chrono::microseconds FunctionTimer::GetAverageGroupTime(const char *szName) {

	GroupDataMap::iterator iter;

	iter = m_mapGroupData.find(szName);

	// Retrieve existing group record
	if (iter != m_mapGroupData.end()) {
		const TimerGroupData & tgd = iter->second;

		return (tgd.iTotalTime / tgd.nEntries);

	// Group record does not exist
	} else {
		return std::chrono::microseconds(0);
	}
}

///////////////////////////////////////////////////////////////////////////////

unsigned long FunctionTimer::GetNumberOfEntries(const char *szName) {

	GroupDataMap::iterator iter;

	iter = m_mapGroupData.find(szName);

	// Retrieve existing group record
	if (iter != m_mapGroupData.end()) {
		const TimerGroupData & tgd = iter->second;

		return (tgd.nEntries);

	// Group record does not exist
	} else {
		return 0;
	}
}

///////////////////////////////////////////////////////////////////////////////

void FunctionTimer::ResetGroupTimeRecord(const char *szName) {
	GroupDataMap::iterator iter;

	iter = m_mapGroupData.find(szName);

	// Retrieve existing group record
	if (iter != m_mapGroupData.end()) {
		m_mapGroupData.erase(iter);

	// Group record does not exist
	} else {
		_EXCEPTION1("Group time record %s does not exist.", szName);
	}
}

///////////////////////////////////////////////////////////////////////////////

