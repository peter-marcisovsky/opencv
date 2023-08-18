#include <iostream>
#include <opencv2/calib3d.hpp>

namespace cv{

int64_t times[50];
std::string tags[50];
int8_t counter = 0;

void saveTime(int64_t time, const char *tag) {
    if (counter < 50) {
        times[counter] = time/1000;
        tags[counter] = std::string(tag);
        counter++;
    }
}

void printTimes() {
    int64_t period = times[counter - 1] - times[0];
    std::cout << "Period = " << period << std::endl << "Intervals" << std::endl;
    for (int i = 0; i < counter - 1; i++) {
        int64_t time_diff = times[i + 1] - times[i];
        std::cout << time_diff << "  \t   " << ((float)time_diff/(float)period)*100 << "% \t\t" << tags[i] << "\t\t <---> \t\t" << tags[i + 1] << std::endl;
    }
}

void clearTimes() {
    for (int i = 0; i < 50; i++) {
        times[i] = 0;
        tags[i].clear();
    }
    counter = 0;
}
}