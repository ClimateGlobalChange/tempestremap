///////////////////////////////////////////////////////////////////////////////
///
///	\file    FunctionTimer.h
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

#ifndef _FUNCTIONTIMER_H_
#define _FUNCTIONTIMER_H_

///////////////////////////////////////////////////////////////////////////////

#include <string>
#include <map>
#include <chrono>

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		FunctionTimer is a class used for timing operations or groups of
///		operations.  Timing is provided via the Unix gettimeofday class
///		and is calculated in microseconds.
///	</summary>
class FunctionTimer {

public:
	///	<summary>
	///		A structure for storing group data.
	///	</summary>
	struct TimerGroupData {
		std::chrono::microseconds iTotalTime;
		unsigned long nEntries;
	};

	///	<summary>
	///		The map structure for storing group data.
	///	</summary>
	typedef std::map<std::string, TimerGroupData> GroupDataMap;

	typedef std::pair<std::string, TimerGroupData> GroupDataPair;

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	///	<param name="szGroup">
	///		Group name used for organizing grouped operations.
	///	</param>
	FunctionTimer(const char *szGroup = NULL);

	///	<summary>
	///		Destructor.
	///	</summary>
	virtual ~FunctionTimer() {
		StopTime();
	}

	///	<summary>
	///		Reset the timer.
	///	</summary>
	void Reset();

	///	<summary>
	///		Return the time elapsed since this timer began.
	///	</summary>
	///	<param name="fDone">
	///		If true stores the elapsed time in the group structure.
	///	</param>
	std::chrono::microseconds Time(bool fDone = false);

	///	<summary>
	///		Return the time elapsed since this timer began and store in
	///		group data.
	///	</summary>
	std::chrono::microseconds StopTime();

public:
	///	<summary>
	///		Retrieve a group data record.
	///	</summary>
	static const TimerGroupData & GetGroupTimeRecord(const char *szName);

	///	<summary>
	///		Retrieve the total time from a group data record.
	///	</summary>
	static std::chrono::microseconds GetTotalGroupTime(const char *szName);

	///	<summary>
	///		Retrieve the average time from a group data record.
	///	</summary>
	static std::chrono::microseconds GetAverageGroupTime(const char *szName);

	///	<summary>
	///		Retrieve the number of entries from a group data record.
	///	</summary>
	static unsigned long GetNumberOfEntries(const char *szName);

	///	<summary>
	///		Reset the group data record.
	///	</summary>
	static void ResetGroupTimeRecord(const char *szName);

private:
	///	<summary>
	///		Group data.
	///	</summary>
	static GroupDataMap m_mapGroupData;

private:
	///	<summary>
	///		Timer is stopped.
	///	</summary>
	bool m_fStopped;

	///	<summary>
	///		Time at which this timer was constructed.
	///	</summary
	std::chrono::high_resolution_clock::time_point m_tpStartTime;

	///	<summary>
	///		Time at which this timer was stopped.
	///	</summary
	std::chrono::high_resolution_clock::time_point m_tpStopTime;

	///	<summary>
	///		Group name associated with this timer.
	///	</summary>
	std::string m_strGroup;
};

///////////////////////////////////////////////////////////////////////////////

#endif

