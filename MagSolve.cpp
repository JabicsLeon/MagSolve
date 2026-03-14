#include "LeoLib.hpp"

#include <array>
#include <exception>
#include <future>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <thread>

//==================================================================================================================Global_struct================================================================================================================
struct Var
{
	long FIELD;
	double var_field;
	double time;
	//long double time_abs;
	std::array<int, 3> date;

	int QMC;
	int ST;

	std::string DATE;
	std::string TIME;		
};

struct var_station 
{
	std::string file_name;

	std::vector<Var> var;
	double X;
	double Y;

	double Time_start;
	double Time_stop;

	/*long double Time_start_abs;
	long double Time_stop_abs;*/

	double dt = 0;
	double mean = 0;

	bool mean_init_ = false;
	bool dt_init_ = false;

	
	std::array<int, 3> Date_start;
	std::array<int, 3> Date_stop;
};

struct Meas
{
	double pk;
	double X;
	double Y;

	double T_top;
	double T_bot;
	double T_grad;
	double T_bot_anom;
	double T_top_anom;


	double dT_var;
	double dT_norm;

	double dT_top;
	double dT_bot;

	double time;
	//long double time_abs;

	int Line;
	int pr;

	std::string DATE;
	std::string TIME;

	std::array<int, 3> date;
};

struct survey 
{
	std::string file_name;

	std::vector<Meas> meas;

	double Time_start;
	double Time_stop;
	
	double mean_bot = 0;
	double mean_top = 0;

	/*long double Time_start_abs;
	long double Time_stop_abs;*/

	std::array<int, 3> Date_start;
	std::array<int, 3> Date_stop;

	bool T_grad_init_ = false;
	bool T_anom_init_ = false;
	bool dT_var_init_ = false;
};

struct configuration
{
	std::vector<survey> sur;
	std::vector<var_station> var_st;

	bool date_init_ = false;
};

//=================================================================================================================Variation_block===============================================================================================================

std::array<size_t, 2> Find_Time(var_station& st, double time)
{
	double dt;
	size_t l = st.var.size(); 

	if ( !st.dt_init_ )
	{
		double t_min = st.var[0].time;
		double t_max = st.var[1].time;
		dt = (t_max - t_min);
		st.dt = dt;
	} 
	else dt = st.dt;
		
	size_t k = (time - st.var[0].time) / dt;

	if ( st.var[k].time == time )
	{
		return { k, k - 1};
	}
	else if ( st.var[k].time >= time )
	{
		while ( k != 0 && ( st.var[k].time >= time && st.var[k - 1].time < time ) != true )
		{
			k -= 1;
		}
	}
	else if ( st.var[k].time <= time )
	{
		while ( k != l - 1 && ( st.var[k].time >= time && st.var[k - 1].time < time ) != true )
		{
			k += 1;
		}
	}

	if ( !( st.var[k].time >= time && st.var[k - 1].time < time ) ) throw std::out_of_range( "Find_Time to variation station: needed time not found!" );

	std::array<size_t, 2> result = { k, k - 1};

	return result;
}


double linear_interpolation(double x3, double x1, double y1, double x2, double y2)
{
	double l = x2 - x1;
	double m = y2 -y1;

	double y3 = m * (x3 - x1) / l + y1;

	return y3;
}


std::vector<double> plane(var_station& st1, var_station& st2, var_station& st3, double time)
{
	std::array<std::array<size_t, 2>, 3> t;
	std::array<var_station, 3> vr = { st1, st2, st3 };

	for (size_t i = 0; i < 3; ++i)
	{
		t[i] = Find_Time(vr[i], time);
	}

	double T1 = linear_interpolation(time, st1.var[t[0][0]].time, st1.var[t[0][0]].var_field, st1.var[t[0][1]].time, st1.var[t[0][1]].var_field);
	double T2 = linear_interpolation(time, st2.var[t[1][0]].time, st2.var[t[1][0]].var_field, st2.var[t[1][1]].time, st2.var[t[1][1]].var_field);
	double T3 = linear_interpolation(time, st3.var[t[2][0]].time, st3.var[t[2][0]].var_field, st3.var[t[2][1]].time, st3.var[t[2][1]].var_field);

	std::vector<double> cof = 
	{ 
		(st2.Y - st1.Y) * (T3 - T1) - (T2 - T1) * (st3.Y - st1.Y),
		(-1) * ((st2.X - st1.X) * (T3 - T1) - (T2 - T1) * (st3.X - st1.X)),
		(st2.X - st1.X) * (st3.Y - st1.Y) - (st3.X - st1.X) * (st2.Y - st1.Y)
	};
	
	cof.push_back((cof[0] * st1.X + cof[1] * st1.Y + cof[2] * T1) * (-1));

	return cof;
	
}


double T_var(survey& sur, var_station& st1, var_station& st2, var_station& st3, size_t k)
{
	double time = sur.meas[k].time;

	std::array<var_station, 3> vr = { st1, st2, st3 };
	std::vector<double> cof = plane(st1, st2, st3, time);

	double T_v = (-1) * ( cof[0] * sur.meas[k].X + cof[1] * sur.meas[k].Y + cof[3] ) / cof[2];

	return T_v;

}


double dT_var(survey& sur, var_station& st1, var_station& st2, var_station& st3, size_t k)
{
	double time = sur.meas[k].time;
	//var_station vr[3] = { st1, st2, st3 };
	std::array<var_station, 3> vr = { st1, st2, st3 };
	std::vector<double> cof = plane(st1, st2, st3, time);
	
	double mean = 0;
	for (size_t i = 0; i < 3; ++i)
	{
		if (!vr[i].mean_init_)
		{
			double tmp = 0;
			for (size_t j = 0; j < vr[i].var.size(); ++j)
			{
				tmp += vr[i].var[j].var_field;
			}
			vr[i].mean_init_ = true;
			vr[i].mean = tmp;
		}

		mean += vr[i].mean;
	}

	mean /= 3;

	double dT_v = (-1) * ( cof[0] * sur.meas[k].X + cof[1] * sur.meas[k].Y + cof[3] ) / cof[2] - mean;

	return dT_v;
	
}


double T_var(survey& sur, var_station& st1, var_station& st2, size_t k)
{
	double time = sur.meas[k].time;
	double l, m, n;

	l = st2.X - st1.X;
	m = st2.Y - st1.Y;

	std::array<std::array<size_t, 2>, 2> t;
	std::array<var_station, 2> vr = { st1, st2 };
	for (size_t i = 0; i < 2; ++i)
	{
		t[i] = Find_Time(vr[i], time);
	}

	double T1 = linear_interpolation(time, st1.var[t[0][0]].time, st1.var[t[0][0]].var_field, st1.var[t[0][1]].time, st1.var[t[0][1]].var_field);
	double T2 = linear_interpolation(time, st1.var[t[1][0]].time, st1.var[t[1][0]].var_field, st1.var[t[1][1]].time, st1.var[t[1][1]].var_field);

	n = T2 - T1;

	leo::matrix<double> S(1,2);

	S[0][0] = l;    S[0][1] = m;

	leo::matrix<double> K(1,2);

	K[0][0] = sur.meas[k].X - st1.X;         K[0][1] = sur.meas[k].Y - st1.Y;

	leo::matrix<double> S_t = S.transposition();

	double T_v = T1 + ( K(S_t) ).get_value() / ( S(S_t) ).get_value();

	return T_v;
}




double dT_var(survey& sur, var_station& st1, var_station& st2, size_t k)
{
	double time = sur.meas[k].time;
	double l, m, n;

	l = st2.X - st1.X;
	m = st2.Y - st1.Y;

	std::array<std::array<size_t, 2>, 2> t;
	std::array<var_station, 2> vr = { st1, st2 };
	for (size_t i = 0; i < 2; ++i)
	{
		t[i] = Find_Time(vr[i], time);
	}


	double T1 = linear_interpolation(time, st1.var[t[0][0]].time, st1.var[t[0][0]].var_field, st1.var[t[0][1]].time, st1.var[t[0][1]].var_field);
	double T2 = linear_interpolation(time, st1.var[t[1][0]].time, st1.var[t[1][0]].var_field, st1.var[t[1][1]].time, st1.var[t[1][1]].var_field);

	n = T2 - T1;

	leo::matrix<double> S(1,2);

	S[0][0] = l;	S[0][1] = m;	//S[0][2] = n;

	leo::matrix<double> K(1,2);

	K[0][0] = sur.meas[k].X - st1.X;	 K[0][1] = sur.meas[k].Y - st1.Y;

	leo::matrix<double> S_t = S.transposition();

	double mean = 0;
	for (size_t i = 0; i < 2; ++i)
	{
		if (!vr[i].mean_init_)
		{
			double tmp = 0;
			for (size_t j = 0; j < vr[i].var.size(); ++j)
			{
				tmp += vr[i].var[j].var_field;
			}
			vr[i].mean_init_ = true;
			vr[i].mean = tmp;
		}
		mean += vr[i].mean;
	}

	mean /= 2;

	double dT_v = T1 + ( K(S_t) ).get_value() / ( S(S_t) ).get_value() - mean;

	return dT_v;
}


double T_var(survey& sur, var_station& st1, size_t k)
{
	double time = sur.meas[k].time;

	std::array<size_t, 2> t = Find_Time(st1, time);

	double T1 = linear_interpolation(time, st1.var[t[0]].time, st1.var[t[0]].var_field, st1.var[t[1]].time, st1.var[t[1]].var_field);

	return T1;
}


double dT_var(survey& sur, var_station& st1, size_t k)
{
	double time = sur.meas[k].time;

	if (!st1.mean_init_)
	{
		double tmp = 0;
		for (size_t j = 0; j < st1.var.size(); ++j)
		{
			tmp += st1.var[j].var_field;
		}
		st1.mean_init_ = true;
		st1.mean = tmp;
	}

	std::array<size_t, 2> t = Find_Time(st1, time);

	double T1 = linear_interpolation(time, st1.var[t[0]].time, st1.var[t[0]].var_field, st1.var[t[1]].time, st1.var[t[1]].var_field);

	double dT_v = T1 - st1.mean;

	return dT_v;
}
//=================================================================================================================Normal_field===============================================================================================================

std::vector<double> T_anom_var(survey& sur, var_station& st1, var_station& st2, var_station& st3, size_t k)
{
	std::vector<double> T_anom;

	double T_vars = T_var(sur, st1, st2, st3, k);

	T_anom.push_back( sur.meas[k].T_bot - T_vars );

	if (sur.T_grad_init_) T_anom.push_back( sur.meas[k].T_top - T_vars );

	return T_anom;
}

std::vector<double> T_anom_var(survey& sur, var_station& st1, var_station& st2, size_t k)
{
	std::vector<double> T_anom;

	double T_vars = T_var(sur, st1, st2, k);

	T_anom.push_back( sur.meas[k].T_bot - T_vars );

	if (sur.T_grad_init_) T_anom.push_back( sur.meas[k].T_top - T_vars );

	return T_anom;
}


std::vector<double> T_anom_var(survey& sur, var_station& st1, size_t k)
{
	std::vector<double> T_anom;

	double T_vars = T_var(sur, st1, k);

	T_anom.push_back( sur.meas[k].T_bot - T_vars );

	if (sur.T_grad_init_) T_anom.push_back( sur.meas[k].T_top - T_vars );

	return T_anom;
}

//=================================================================================================================Leveling_data===============================================================================================================
void Leveling(configuration config)
{
	double mean_bot = 0;
	double mean_top = 0;

	for(auto& it : config.sur)
	{
		mean_bot += it.mean_bot;
		if (it.T_grad_init_) mean_top += it.mean_top;
	}

	double N = config.sur.size();
	mean_bot /= N;
	mean_top /= N;

	std::vector<std::thread> proc;

	for(auto& it : config.sur)
	{
		proc.push_back( std::thread( [&it, mean_bot, mean_top]
		{
			double dmean_bot = mean_bot - it.mean_bot;
			double dmean_top = 0;

			if (it.T_grad_init_)
			{
				 dmean_top = mean_top - it.mean_top;
			}

			if (it.T_anom_init_)
			{
				if (it.T_grad_init_)
				{
					for(auto& el : it.meas)
					{
						el.T_bot_anom += dmean_bot;
						el.T_top_anom += dmean_top;
					}
				}
				else
				{
					for(auto& el : it.meas)
					{
						el.T_bot_anom += dmean_bot;
					}
				}
			}
		}
		));
	}

	for(auto& t : proc)
	{
		t.join();
	}
}

//=================================================================================================================Parsing_File===============================================================================================================


var_station VarParser(std::stringstream& rd, const std::string& delimiter="\t")
{
	var_station vr;
	double mean = 0;
	double dt = 0;

	std::string line;
	size_t i = 0;
	std::getline(rd, line);
	while (std::getline(rd, line))
	{
		std::vector<std::string> date = leo::split(line, delimiter);

		if (date.size() < 5)
		{
			std::cerr << "Mistake in line: " << i << "\n";
			throw std::out_of_range("VarParser: There is string in ypu file, that has less number of argument them excepted (5)");
		}

		vr.var.resize(i + 1);

		vr.var[i].FIELD = leo::string_to<long>(date[0]);
		vr.var[i].QMC = leo::string_to<int>(date[1]);
		vr.var[i].ST = leo::string_to<int>(date[2]);
		vr.var[i].DATE = date[3];
		vr.var[i].TIME = date[4];

		vr.var[i].var_field = double(vr.var[i].FIELD) / 1000;
		mean += vr.var[i].var_field;

		std::vector<std::string> tm = leo::split(vr.var[i].TIME, ":");
		if (tm.size() !=3) 
		{
			std::cerr << "Mistake in line: " << i + 1<< "\n";
			std::cerr << "Invalid value: " << vr.var[i].TIME << "\n";
			throw std::invalid_argument("VarParser: Invalid format of TIME");
		}
		vr.var[i].time = leo::string_to<double>(tm[0]) * 60 * 60 +
				leo::string_to<double>(tm[1]) * 60 +
				leo::string_to<double>(tm[2]);

		std::vector<std::string> date_tm = leo::split(vr.var[i].DATE, ".");
		if (date_tm.size() !=3)
		{
			std::cerr << "Mistake in line: " << i << "\n";
			throw std::invalid_argument("VarParser: Invalid format of DATE");
		}
		vr.var[i].date[0] = leo::string_to<int>(date_tm[0]);
		vr.var[i].date[1] = leo::string_to<int>(date_tm[1]);
		vr.var[i].date[2] = leo::string_to<int>(date_tm[2]);

		++i;
	}

	vr.mean = mean / i;
	vr.dt = ( vr.var[1].time - vr.var[0].time );
	vr.mean_init_ = true;
	vr.dt_init_ = true;
	
	vr.Time_start = vr.var[0].time;
	vr.Time_stop = vr.var[i - 1].time;

	vr.Date_start = vr.var[0].date;
	vr.Date_stop = vr.var[i - 1].date;	

	return vr;	
}

survey SurParser(std::stringstream& rd, const std::string& delimiter="\t")
{
	survey sur;
	
	std::string line;
	size_t i = 0;
	std::getline(rd, line);
	while (std::getline(rd, line))
	{
		std::vector<std::string> date = leo::split(line, delimiter);

		if (date.size() < 8)
		{
			std::cerr << "Mistake in line: " << i << "\n";
			std::cerr << "Your line: " << line << "\n";
			throw std::out_of_range("SurParser: There is string in your file, that has less number of argument them excepted (9)");

		}

		sur.meas.resize(i + 1);
		
		sur.meas[i].Line = leo::string_to<int>(date[0]);
		sur.meas[i].pr = leo::string_to<int>(date[1]);
		sur.meas[i].pk = leo::string_to<double>(date[2]);
		sur.meas[i].X = leo::string_to<double>(date[3]);
		sur.meas[i].Y = leo::string_to<double>(date[4]);
		sur.meas[i].TIME = date[5];
		sur.meas[i].DATE = date[6];
		sur.meas[i].T_bot = leo::string_to<double>(date[7]) / 1000;
		sur.mean_bot += sur.meas[i].T_bot;
		if (date.size() > 8)
		{
			sur.meas[i].T_top = leo::string_to<double>(date[8]) / 1000;
			sur.meas[i].T_grad = sur.meas[i].T_top - sur.meas[i].T_bot;
			sur.mean_top += sur.meas[i].T_top;
			if (!sur.T_grad_init_) sur.T_grad_init_ = true;
		}

		std::vector<std::string> tm = leo::split(sur.meas[i].TIME, ":");
		if (tm.size() !=3)
		{
			std::cerr << "Mistake in line: " << i << "\n";
			std::cerr << "Invalid value: " << sur.meas[i].TIME << "\n";
			throw std::invalid_argument("SurParser: Invalid format of TIME");
		}
		sur.meas[i].time = leo::string_to<double>(tm[0]) * 60 * 60 +
				leo::string_to<double>(tm[1]) * 60 +
				leo::string_to<double>(tm[2]);

		std::vector<std::string> date_tm = leo::split(sur.meas[i].DATE, ".");
		if (date_tm.size() !=3)	
		{
			std::cerr << "Mistake in line: " << i << "\n";
			throw std::invalid_argument("VarParser: Invalid format of DATE");
		}
		sur.meas[i].date[0] = leo::string_to<int>(date_tm[0]);
		sur.meas[i].date[1] = leo::string_to<int>(date_tm[1]);
		sur.meas[i].date[2] = leo::string_to<int>(date_tm[2]);

		++i;
	}

	sur.mean_bot /= sur.meas.size();
	if(sur.T_grad_init_) sur.mean_top /= sur.meas.size();

	sur.Time_start = sur.meas[0].time;
	sur.Time_stop = sur.meas[i - 1].time;

	sur.Date_start = sur.meas[0].date;
	sur.Date_stop = sur.meas[i - 1].date;

	return sur;
}

std::stringstream SurWrite(survey sur, const std::string& del="\t", bool head=true)
{
	std::stringstream ss;
	ss << std::fixed << std::setprecision(4);

	if (head)
	{
		ss << "Line" << del
			 << "pr" << del
			 << "pk" << del
			 << "X_end" << del
			 << "Y_end" << del
			 << "TIME" << del
			 << "DATE";

		if (sur.T_grad_init_) ss << del 
					<< "T_bottom" << del
					<< "T_top" << del
					<< "T_grad";
		else ss << del << "T_measurment";

		if (sur.dT_var_init_) 
		{
			ss << del << "dT_variation";

			if (sur.T_grad_init_) ss << del 
						<< "dT_bottom" << del
						<< "dT_top";
			else ss << del << "dT";
		}

		if (sur.T_anom_init_)
		{
			if (sur.T_grad_init_) ss << del
						<< "T_bottom_anomal" << del
						<< "T_top_anomal";
			else ss << "T_anomal";
		}

		ss << "\n";
	}


	for (auto& it : sur.meas)
	{
		ss << it.Line << del
			<< it.pr << del
			<< it.pk << del
			<< it.X << del
			<< it.Y << del
			<< it.TIME << del
			<< it.DATE << del
			<< it.T_bot;

		if (sur.T_grad_init_) ss << del
						<< it.T_top << del
						<< it.T_grad;

		if (sur.dT_var_init_) 
		{
			ss << del << it.dT_var;

			if (sur.T_grad_init_) ss << del 
						<< it.dT_bot << del
						<< it.dT_top;
			else ss << del << it.dT_bot;
		}

		if (sur.T_anom_init_)
		{
			if (sur.T_grad_init_) ss << del
						<< it.T_bot_anom << del
						<< it.T_top_anom;
			else ss << it.T_bot_anom;
		}

		ss << "\n";
	}

	return ss;
}


configuration ConfParser(std::stringstream& rd, const std::string& delimiter="\t")
{
	configuration config;

	std::string line;
	size_t i = 0;
	while (std::getline(rd, line))
	{
		std::vector<std::string> date = leo::split(line, delimiter);
		
		if ( date[0] == "-var" || date[0] == "--variation" )
		{
			if (date.size() < 5)
			{
				std::cerr << "Error of configuration file in line " << i + 1 << ": excepted 5 argument in line!\n";
				std::cerr << "Correct forman: [-var/--variation]   [[-n/--name] [file_name]]   [-x=.../-X=...]   [-y=.../-Y=...]\n";
				throw std::invalid_argument("ConfParser: Error format of file");
			}

			std::string file_name;
			double X, Y;
			bool X_init_ = false;
			bool Y_init_ = false;
			bool file_name_init_ = false;

			size_t j = 1;
			while ( j < date.size() )
			{
				if ( date[j] == "-n" || date[j] == "--name" )
				{
					file_name = date[j + 1];
					file_name_init_ = true;
					++j;
				} 
				else 
				{
					std::vector<std::string> cord = leo::split(date[j], "=");
					if ( cord.size() < 3 )
					{
						if ( cord[0] == "-x" || cord[0] == "-X" )
						{
							X = leo::string_to<double>(cord[1]);
							X_init_ = true;
						}
						else if ( cord[0] == "-y" || cord[0] == "-Y" )
						{
							Y = leo::string_to<double>(cord[1]);
							Y_init_ = true;
						}
						else
						{
							std::cerr << "Error of configuration file in line " << i + 1 << ": unexcepted value: [" << date[j] << "]\n";
							std::cerr << "Correct forman to variation station: [-var/--variation]   [[-n/--name] [file_name]]   [-x=.../-X=...]   [-y=.../-Y=...]\n";
							throw std::invalid_argument("ConfParser: Error format of file");
						}
					}
				}
				++j;
			}
			
			if ( !X_init_ || !Y_init_ || file_name_init_ )
			{
				if (!X_init_)
				{
					std::cerr << "Error of configuration file in line " << i + 1 << ": X was't init!\n";
					throw std::invalid_argument("ConfParser: Error format of file");
				}
				if (!Y_init_)
				{
					std::cerr << "Error of configuration file in line " << i + 1 << ": Y was't init!\n";
					throw std::invalid_argument("ConfParser: Error format of file");
				}
				if (!file_name_init_)
				{
					std::cerr << "Error of configuration file in line " << i + 1 << ": file way was't init!\n";
				}
			}

			config.var_st.resize(config.var_st.size() + 1);
			config.var_st[config.var_st.size() - 1].file_name = file_name;
			config.var_st[config.var_st.size() - 1].X = X;
			config.var_st[config.var_st.size() - 1].Y = Y;

		}
		else if ( date[0] == "-meas" || date[0] == "--measurment" )
		{
			if (date.size() < 3)
			{
				std::cerr << "Error of configuration file in line " << i + 1 << ": excepted 3 argument in line!\n";
				std::cerr << "Correct forman: [-meas/--measurment] [-n/--name] [file_name]\n";
				throw std::invalid_argument("ConfParser: Error format of file");
			} 
			if ( date[1] == "-n" || date[1] == "--name" )
			{
				config.sur.resize(config.sur.size() + 1);
				config.sur[config.sur.size() - 1].file_name = date[2];

			}
		} 
		else
		{
			std::cerr << "Error of configuration file in line " << i + 1 << ": excepted definition type of file!\n";
			std::cerr << "In your file: [" << date[0] << "]\n";
			std::cerr << "Write [-var] or [--variation] to add variation station\n";
			std::cerr << "Write [-meas] or [--measurment] to add observation\n";
			std::cerr << "Correct forman to variation station: [-var/--variation]   [[-n/--name] [file_name]]   [-x=.../-X=...]   [-y=.../-Y=...]\nCorrect format to observation: [-meas/--measurment]   [[-n/--name] [file_name]]\n";
			throw std::invalid_argument("ConfParser: Error format of file");
		}

		++i;
	}

	return config;
}

//=================================================================================================================Init===============================================================================================================


void ConfigInit(configuration& config, const std::string& delimiter="\t")
{
	size_t len_sur = config.sur.size();
	size_t len_var_st = config.var_st.size();

	std::vector<std::thread> sur_tr;
	std::vector<std::thread> var_st_tr;
	/*std::vectot<std::promis<void>> prom_sur;
	std::vectot<std::promis<void>> prom_var;
	std::vector<std::future<void>> fut_sur;
	std::vector<std::future<void>> fut_var;*/

	for (size_t i = 0; i < len_sur; ++i)
	{
		//prom_sur.resize(prom_sur.size() + 1);
		//fut_var.resize(fut_var.size() + 1);

		//fut_sur[i] = prom_sur[i].get_future();

		sur_tr.push_back( std::thread(  [&config, i, &delimiter/*, &fut_sur[i]*/] 
		{
			std::stringstream ss = leo::ReadFile(config.sur[i].file_name);
			survey sur = SurParser(ss, delimiter);
			sur.file_name = config.sur[i].file_name;
			config.sur[i] = sur;
			//fut_var[i].set_value();
		}
		));
	}

	for (size_t i = 0; i < len_var_st; ++i)
	{
		var_st_tr.push_back( std::thread( [&config, i, &delimiter]
		{
			std::stringstream ss = leo::ReadFile(config.var_st[i].file_name);
			var_station var = VarParser(ss, delimiter);
			var.file_name = config.var_st[i].file_name;
			var.X = config.var_st[i].X;
			var.Y = config.var_st[i].Y;
			config.var_st[i] = var;
		}
		));
	}

	for (auto& t : sur_tr)
	{
		t.join();
	}

	for (auto& t : var_st_tr)
	{
		t.join();
	}

	config.date_init_ = true;
}

double abs_time(const std::array<int, 3>& d, double t)
{
	return (d[2] * 365 + d[1] * 30 + d[0]) * 24 * 60 * 60 + t;
}

bool covers(const std::array<int, 3>& date_start_st, double time_start_st,
		const std::array<int, 3>& date_stop_st, double time_stop_st,
		const std::array<int, 3>& date_start_sur, double time_start_sur,
		const std::array<int, 3>& date_stop_sur, double time_stop_sur)
{
	auto to_abs = [] (const std::array<int, 3>& d, double t)
	{
		return (d[2] * 365 + d[1] * 30 + d[0]) * 24 * 60 * 60 + t;
	};

	double st_start = to_abs(date_start_st, time_start_st);
	double st_stop = to_abs(date_stop_st, time_stop_st);
	double sur_start = to_abs(date_start_sur, time_start_sur);
	double sur_stop = to_abs(date_stop_sur, time_stop_sur);

	return st_start <= sur_start && sur_stop <= st_stop;

}



void dT_varInit(configuration& config)
{
	size_t len = config.sur.size();

	std::vector<std::thread> sur_tr;
	
	for (size_t i = 0; i < len; ++i)
	{
		sur_tr.push_back( std::thread(  [&config, i]
		{
			int j = 0;
			
			std::vector<var_station> vr_st;

			for (auto& it : config.var_st)
			{
				if ( covers(it.Date_start, it.Time_start, it.Date_stop, it.Time_stop,
						config.sur[i].Date_start, config.sur[i].Time_start,
						config.sur[i].Date_stop, config.sur[i].Time_stop) )
						
				{
					vr_st.push_back(it);
					++j;
				}
			
			}


			if (j >= 3)
			{
				for (size_t k = 0; k < config.sur[i].meas.size(); ++k)
				{
					double dT_v = dT_var(config.sur[i], vr_st[0], vr_st[1], vr_st[2], k);
					config.sur[i].meas[k].dT_var = dT_v;

					config.sur[i].meas[k].dT_bot = config.sur[i].meas[k].T_bot - dT_v;

					if (config.sur[i].T_grad_init_)
					{
						config.sur[i].meas[k].dT_top = config.sur[i].meas[k].T_top - dT_v;
					}
				}

				config.sur[i].dT_var_init_ = true;
			}
			else if (j == 2)
			{
				for (size_t k = 0; k < config.sur[i].meas.size(); ++k)
				{
					double dT_v = dT_var(config.sur[i], vr_st[0], vr_st[1], k);
					config.sur[i].meas[k].dT_var = dT_v;

					config.sur[i].meas[k].dT_bot = config.sur[i].meas[k].T_bot - dT_v;

					if (config.sur[i].T_grad_init_)
					{
						config.sur[i].meas[k].dT_top = config.sur[i].meas[k].T_top - dT_v;
					}
				}

				config.sur[i].dT_var_init_ = true;
			}
			else if (j == 1)
			{
				for (size_t k = 0; k < config.sur[i].meas.size(); ++k)
				{
					double dT_v = dT_var(config.sur[i], vr_st[0], k);

					config.sur[i].meas[k].dT_var = dT_v;

					config.sur[i].meas[k].dT_bot = config.sur[i].meas[k].T_bot - dT_v;

					if (config.sur[i].T_grad_init_)
					{
						config.sur[i].meas[k].dT_top = config.sur[i].meas[k].T_top - dT_v;
					}
				}

				config.sur[i].dT_var_init_ = true;
			}
			else
			{
				std::cerr << "dT_varInit: Error of searching station in range time! file:"
						<< config.sur[i].file_name << "\n";
				std::cerr << "Excepted Date range: " << config.sur[i].Date_start[0] << "."
								<< config.sur[i].Date_start[1] << "."
								<< config.sur[i].Date_start[2] << " - "
								<< config.sur[i].Date_stop[0] << "." 
								<< config.sur[i].Date_stop[1] << "."
								<< config.sur[i].Date_stop[2] << "\n"; 
				throw std::invalid_argument("dT_varInit: There is't variation staition in observation time range!");
			}
		}
		));
	}

	for (auto& t : sur_tr)
	{
		t.join();
	}

}



void T_anom_varInit(configuration& config)
{
	size_t len = config.sur.size();

	std::vector<std::thread> sur_tr;

	for (size_t i = 0; i < len; ++i)
	{
		sur_tr.push_back( std::thread(  [&config, i]
		{
			int j = 0;

			std::vector<var_station> vr_st;

			for (auto& it : config.var_st)
			{
				if ( covers(it.Date_start, it.Time_start, it.Date_stop, it.Time_stop,
						config.sur[i].Date_start, config.sur[i].Time_start,
						config.sur[i].Date_stop, config.sur[i].Time_stop) )
				{
					vr_st.push_back(it);
					++j;
				}
			}

			if (j >= 3)
			{
				for (size_t k = 0; k < config.sur[i].meas.size(); ++k)
				{
					std::vector<double> T_a = T_anom_var(config.sur[i], vr_st[0], vr_st[1], vr_st[2], k);

					config.sur[i].meas[k].T_bot_anom = T_a[0];
					if (config.sur[i].T_grad_init_) config.sur[i].meas[k].T_top_anom = T_a[1];
				}

				config.sur[i].T_anom_init_ = true;
			} 
			else if (j == 2) 
			{
				for (size_t k = 0; k < config.sur[i].meas.size(); ++k)
				{
					std::vector<double> T_a = T_anom_var(config.sur[i], vr_st[0], vr_st[1], k);

					config.sur[i].meas[k].T_bot_anom = T_a[0];
					if (config.sur[i].T_grad_init_) config.sur[i].meas[k].T_top_anom = T_a[1];
				}

				config.sur[i].T_anom_init_ = true;
			}
			else if (j == 1)
			{
				for (size_t k = 0; k < config.sur[i].meas.size(); ++k)
				{
					std::vector<double> T_a = T_anom_var(config.sur[i], vr_st[0], k);

					config.sur[i].meas[k].T_bot_anom = T_a[0];
					if (config.sur[i].T_grad_init_) config.sur[i].meas[k].T_top_anom = T_a[1];
				}

				config.sur[i].T_anom_init_ = true;
			}
			else
			{
				std::cerr << "dT_varInit: Error of searching station in range time! file:"
						 << config.sur[i].file_name << "\n";
				std::cerr << "Excepted Date range: " << config.sur[i].Date_start[0] << "."
						<< config.sur[i].Date_start[1] << "."
						<< config.sur[i].Date_start[2] << " - "
						<< config.sur[i].Date_stop[0] << "."
						<< config.sur[i].Date_stop[1] << "."
						<< config.sur[i].Date_stop[2] << "\n";
				throw std::invalid_argument("dT_varInit: There is't variation staition in observation time range!");
			}
		}
		));
	}

	for (auto& t : sur_tr)
	{
		t.join();
	}
}



std::stringstream VarWrite(var_station& st, const std::string& del="\t")
{
	std::stringstream ss;
	ss << std::fixed << std::setprecision(4);

	ss << "FIELD" << del
		<< "QMC" << del
		<< "ST" << del
		<< "DATE" << del
		<< "TIME" << del
		<< "var_field" << del 
		<< "time [sec]" << del
		<< "abs time [sec]" << "\n";

	//ss << std::fixed << std::setprecision(4);

	for (auto& it : st.var)
	{
		ss << it.FIELD << del
			<< it.QMC << del
			<< it.ST << del
			<< it.DATE << del
			<< it.TIME << del
			<< it.var_field << del
			<< it.time << del
			<< abs_time(it.date, it.time) << "\n";
	}

	return ss;
}





void CorrectFormatInput()
{
	std::cerr << "Correct format:\n"
		<< "[comands...] [-config/--configuration] [config_file] [options...]\n"
		<< "Comands:\n"
		<< "\t[-var/--variation] - calculate variation by station\n"
		<< "\t[-pvar/--printVar] - calculate and print variation station whit time in second\n"
		<< "\t[-avar/--anomVar] - calculate anomals field by field of variation station\n"
		<< "\t[-lev/--leveling] - leveling all date\n"
		<< "Options:\n"
		<< "\t[-del/--delimiter] [symbol] - set delimiter in reads files\n"
		<< "\t[-of/--outfile] - creat output files whit prefix 'processing_'\n"
		<< "\t[-of/--outfile] [file] - write output to file\n"
		<< "\t[-of/--outfile] [file] [-un/--unitMeas] - join all observation";
}

void CorrectFormatConfigFile()
{
	std::cerr << "Correct format of line in [config_file]:\n"
		<< "Forman to variation station:\n"
		<< "\t[-var/--variation] [-n/--name] [file_name] [-x=.../-X=...]   [-y=.../-Y=...]\n"
		<< "\tSystem UTM\n"
		<< "Format to observation:\n"
		<< "\t [-meas/--measurment] [-n/--name] [file_name]\n";
}


std::string basename(const std::string& path)
{
	size_t pos = path.find_last_of("/\\");
	if (pos == std::string::npos) return path;
	return path.substr(pos + 1);
}

const std::vector<std::string> names =
{
	"-var",
	"--variation",
	"-config",
	"--configuration",
	"-of",
	"--outfile",
	"-del",
	"--delimiter",
	"-pvar",
	"--printVar",
	"-avar",
	"--anomVar",
	"-un",
	"--unitMeas",
	"-lev",
	"--leveling"
};

bool not_option(std::string com)
{
	bool result = true;

	for (const auto& it : names)
	{
		result = result && (com != it);
	}

	return result;
}


//=================================================================================================================Parsing_Dialog===============================================================================================================



int main(int argc, char* argv[]) 
{
	bool vars = false;
	bool to_file = false;
	bool create_file = false;
	bool print_var = false;
	bool T_an_var = false;
	bool Unit_meas = false;
	bool Level = false;

	std::string delimiter = "\t";
	std::string config_file;
	std::string out_file;

	if (argc < 4)
	{
		std::cerr << "Invalid input format! Excepted number of argument: more then 3. You have inputed: " << argc - 1 << "\n";
		CorrectFormatInput();
		std::cerr << "\n";
		CorrectFormatConfigFile();
		return 1;
	}

	
	for (size_t i = 1; i < argc; ++i)
	{
		std::string arg1 = argv[i];
		std::string arg2 = "";
		if (argc - i >= 2) arg2 = argv[i + 1];
		if (arg1 == "-var" || arg1 == "--variation")
		{
			vars = true;
		}
		else if (arg1 == "-config" || arg1 == "--configuration")
		{
			if (argc - i >= 2)
			{
				config_file = arg2;
				++i;
			}
			else
			{
				std::cerr << "Invalid input format!";
				CorrectFormatInput();
			}
		}
		else if (arg1 == "-of" || arg1 == "--outfile")
		{
			if (argc - i >= 2 && not_option(arg2))
			{
				out_file = argv[i + 1];
				to_file = true;
				create_file = true;
				++i;
			}
			else
			{
				create_file = true;
			}
		}
		else if (arg1 == "-del" || arg1 == "--delimiter")
		{
			if (argc - i >= 2)
			{
				delimiter = arg2;
			}
			else
			{
				std::cerr << "Invalid input format!";
				CorrectFormatInput();
			}
		}
		else if (arg1 == "-pvar" || arg1 == "--printVar")
		{
			print_var = true;
		}
		else if (arg1 == "-avar" || arg1 == "--anomVar")
		{
			T_an_var = true;
		}
		else if (arg1 == "-un" || arg1 == "--unitMeas")
		{
			Unit_meas = true;
		}
		else if (arg1 == "-lev" || arg1 == "--leveling")
		{
			Level = true;
		}
		else
		{
			std::cerr << "Invalid input format! Unexcepted argument: " << argv[i] << "\n";
			CorrectFormatInput();
			return 1;
		}
	}


	std::stringstream confss = leo::ReadFile(config_file);
	configuration config = ConfParser(confss, delimiter);
	ConfigInit(config, delimiter);

	std::vector<std::thread> processing;

	if (vars) processing.push_back( std::thread(dT_varInit, std::ref(config)) );
	if (T_an_var) processing.push_back( std::thread(T_anom_varInit, std::ref(config)) );

	for(auto& t : processing)
	{
		t.join();
	}


	if (Level) Leveling(config);

	
	if (create_file)
	{
		if (to_file)
		{
			std::stringstream ss;
			bool _init_head_ = false;			

			for (auto& it : config.sur)
			{
				bool phead = Unit_meas && !(_init_head_);
				_init_head_ = true;

				if (!Unit_meas) ss << it.file_name << "\n";
				ss << std::fixed << std::setprecision(4)
					<< SurWrite(it, "\t", phead).str();
				ss << "\n\n";
			}
	
			if (print_var)
			{
				for (auto& it : config.var_st)
				{
					ss << it.file_name << "\n";
					ss << std::fixed << std::setprecision(4)
						<< VarWrite(it).str();
					ss << "\n\n";
				}
			}


			leo::WriteFile(out_file, ss);
		}
		else
		{
			for (auto& it : config.sur)
			{
				std::stringstream ss;

				std::string file_name = "processing_" + basename(it.file_name);
				
				ss << std::fixed << std::setprecision(4)
					<< SurWrite(it).str();

				leo::WriteFile(file_name, ss);
			}


			if (print_var)
			{
				for (auto& it : config.var_st)
				{
					std::stringstream ss;

					std::string file_name = "processing_" + basename(it.file_name);
					
					ss << std::fixed << std::setprecision(4)
						<< VarWrite(it).str();

					leo::WriteFile(file_name, ss);
				}
			}
		}
	}
	else
	{
		std::stringstream ss;
		bool _init_head_ = false;

		for (auto& it : config.sur)
		{
			bool phead = Unit_meas && !(_init_head_);
			_init_head_ = true;

			if (!Unit_meas) ss << it.file_name << "\n";
			ss << std::fixed << std::setprecision(4)
				<< SurWrite(it, "|", phead).str();
			ss << "\n\n";
		}


		if (print_var)
		{
			for (auto& it : config.var_st)
			{
				ss << it.file_name << "\n";
				ss << std::fixed << std::setprecision(4)
					<< VarWrite(it, "|").str();
				ss << "\n\n";
			}
		}


		std::cout << std::fixed << std::setprecision(4)
			<< ss.str();
	}

	std::cout << "Creat by LeoLib. MSU 2026\n";

	return 0;
}
