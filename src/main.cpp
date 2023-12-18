#include <sim.hpp>
#include <nlohmann/json.hpp>

int main(int argc, char** argv)
{
    // Set hyper params
    ecl::device_count = 200;
    ecl::file_count = 100;
    ecl::epsilon = 1e-4;
    ecl::cs = 1;
    ecl::Q = 20;
    ecl::S = 50;
    ecl::c0 = 0.3;
    double demand_factor = 30;
    double discount_factor = 0.1;
    double contant_sharing_price = 0.1;
    double update_rate = 0.05;
    ecl::zipf zipf_dist(ecl::file_count);
    double lambda = 0.4;

    // Create devices and files
    std::vector<ecl::Device> devices;
    std::vector<ecl::File> files;
    
    for(int i = 0; i < ecl::device_count; i++)
        devices.push_back(ecl::Device(demand_factor * ecl::uniform(ecl::gen), ecl::uniform(ecl::gen), ecl::Q, ecl::S, i, ecl::file_count));

    for(int i = 0; i < ecl::file_count; i++)
        files.push_back(ecl::File(0, zipf_dist(i+1), contant_sharing_price+discount_factor*zipf_dist(i+1)));


    // Compute init fi and alpha
    
    std::vector<double> fis;
    std::vector<double> alphas;
    fis.resize(ecl::file_count);
    alphas.resize(ecl::file_count);

    // while loop
    int test = 0;
    printf("\n");
    auto bar = [](int a, int b)->std::string{
        int barLength = 50;
        int stop = barLength * (static_cast<double>(a) / static_cast<double>(b));
        std::string s;
        for(int i = 0; i < stop; i++) s.append("=");
        for(int i = 0; i < barLength-stop; i++) s.append(".");
        return s;
    };
    
    // logging
    std::fstream loggingfile;
    auto now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    std::string logLocation;//std::put_time(std::localtime(&now), "%Y.%m.%d_%H:%M:%S");
    std::stringstream tts;
    tts << std::put_time(std::localtime(&now), "%Y.%m.%d_%H:%M:%S");
    logLocation = "logs/" + tts.str() + ".log";
    loggingfile.open(logLocation, std::ios::out);
    if(!loggingfile.is_open())
        exit(EXIT_FAILURE);
    std::vector<int> f1s, f2s, f3s, f4s;    

    int epoch = 200;
    int f1=0, f2=0, f3=0, f4=0;
    while(test++<epoch)
    {
        printf("\033[1AEpoch:%d/%d <%s>\n",test, epoch, bar(test, epoch).c_str());
        std::vector<std::vector<int>> strategies;
        strategies.resize(ecl::device_count);
        for(int i = 0; i < ecl::device_count; i++)
            strategies[i] = devices[i].MEC_PEC_Strategy_Update(alphas, fis, files, devices);
        for(int i = 0; i < ecl::device_count; i++)
            devices[i].caching_desicions = strategies[i];
        
        auto newfis = fis;
        auto newalphas = alphas;
        for(int i = 0; i < ecl::file_count; i++)
        {
            newalphas[i] = ecl::get_alpha(i, ecl::Q, files, devices);
            newfis[i] = lambda * fis[i] + (1 - lambda) * ecl::Dk(i, files, devices) * files[i].price / ecl::iterativeS1(i, files, devices, alphas, fis);//ecl::get_fi(i, ecl::Q, files, devices);
        }
        for(int i = 0; i < ecl::file_count; i++)
        {
            alphas[i] = update_rate * newalphas[i] + (1 - update_rate) * alphas[i];
            fis[i] = update_rate * newfis[i] + (1 - update_rate) * fis[i];
        }

        // statistics
        f1=0, f2=0, f3=0, f4=0;
        for(int i = 0; i < ecl::device_count; i++)
        {
            f1 += (devices[i].caching_desicions[0] == 1);
            f2 += (devices[i].caching_desicions[1] == 1);
            f3 += (devices[i].caching_desicions[2] == 1);
            f4 += (devices[i].caching_desicions[3] == 1);
        }
        f1s.push_back(f1);
        f2s.push_back(f2);
        f3s.push_back(f3);
        f4s.push_back(f4);


    }
    
    nlohmann::json target;
    target["device_count"] = ecl::device_count;
    target["file_count"] = ecl::file_count;
    target["cs"] = ecl::cs;
    target["Q"] = ecl::Q;
    target["S"] = ecl::S;
    target["c0"] = ecl::c0;
    target["demand_factor"] = demand_factor;
    target["discount_factor"] = discount_factor;
    target["constant_sharing_price"] = contant_sharing_price;
    target["epoch"] = epoch;
    target["f1"] = f1s;
    target["f2"] = f2s;
    target["f3"] = f3s;
    target["f4"] = f4s;
    loggingfile << target.dump();
    loggingfile.close();

    f1=0, f2=0, f3=0, f4=0;
    for(int i = 0; i < ecl::device_count; i++)
    {
        f1 += (devices[i].caching_desicions[0] == 1);
        f2 += (devices[i].caching_desicions[1] == 1);
        f3 += (devices[i].caching_desicions[2] == 1);
        f4 += (devices[i].caching_desicions[3] == 1);
    }

    printf("File 1 Agent:%d/%d, %.2lf%%\nFile 2 Agent:%d/%d, %.2lf%%\nFile 3 Agent:%d/%d, %.2lf%%\nFile 4 Agent:%d/%d, %.2lf%%\n",
        f1, ecl::device_count, 100 * f1 / static_cast<double>(ecl::device_count),
        f2, ecl::device_count, 100 * f2 / static_cast<double>(ecl::device_count),
        f3, ecl::device_count, 100 * f3 / static_cast<double>(ecl::device_count),
        f4, ecl::device_count, 100 * f4 / static_cast<double>(ecl::device_count));

    auto t = ecl::tracker;
    return 0;
}