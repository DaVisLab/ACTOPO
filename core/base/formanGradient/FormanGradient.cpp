#include <FormanGradient.h>

FormanGradient::FormanGradient() {
  // logFile_.open("log_fg.txt", std::ofstream::trunc);
  // if(logFile_.fail()) {
  //   this->printErr("Cannot open file for logging!");
  // }
  relationTime_ = 0.0;
  this->setDebugMsgPrefix("FormanGradient");
}

FormanGradient::~FormanGradient() {
  // printMsg("Time spent for relation computation is " + std::to_string(relationTime_));
}

void FormanGradient::setPair(const Simplex &tail, const Simplex &head) {

  if(tail.first == 0) {
    gradient_[tail.first][tail.second][0] = head.second;
  } else {
    gradient_[tail.first][tail.second][1] = head.second;
  }

  gradient_[head.first][head.second][0] = tail.second;
}

bool FormanGradient::getPair(const Simplex &simpl, Simplex &paired) {
  if(simpl.first == 0) {
    if(gradient_[simpl.first][simpl.second][0] != -1) {
      paired = Simplex(1, gradient_[simpl.first][simpl.second][0]);
      return true;
    }
  } else if(simpl.first == dimensionality_) {
    if(gradient_[simpl.first][simpl.second][0] != -1) {
      paired
        = Simplex(dimensionality_ - 1, gradient_[simpl.first][simpl.second][0]);
      return true;
    }
  } else {
    if(gradient_[simpl.first][simpl.second][0] != -1) {
      paired
        = Simplex(simpl.first - 1, gradient_[simpl.first][simpl.second][0]);
      return true;
    } else if(gradient_[simpl.first][simpl.second][1] != -1) {
      paired
        = Simplex(simpl.first + 1, gradient_[simpl.first][simpl.second][1]);
      return true;
    }
  }

  return false;
}

bool FormanGradient::isCritical(const Simplex &simpl) {
  for(auto i : gradient_[simpl.first][simpl.second]) {
    if(i != -1)
      return false;
  }

  return true;
}