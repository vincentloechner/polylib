pipeline {
  agent none
  stages{
    stage('PolyLib'){
      //Using a matrix section for extensibility
      matrix{
        //The work is done on a slave with the indicated OS
        agent{ label "${OS}" }
        axes{
          axis{
            //The list of target OS
            name 'OS'
            values 'Ubuntu', 'macOS', 'CentOS', 'fedora', 'Debian'
          }
        }
        //Here is the actual work to be done
        stages{
          //Install whatever dependency the project may have
          stage('Tools'){
            steps{
              script{
                if(env.OS == 'macOS')
                  sh 'brew install automake libtool'
                if(env.OS == 'CentOS')
                  sh 'sudo yum install gmp-devel -y'
                if(env.OS == 'fedora')
                  sh 'sudo dnf install gmp-devel -y'
                if(env.OS == 'Debian')
                  sh 'sudo apt install autoconf libtool libgmp-dev make -y'
              }
            }
          }
          //Build the project
          stage('Build'){
            steps{
              sh './autogen.sh && ./configure && make -j'
            }
          }
          //Execute the test suites
          stage('Test'){
            steps{
              sh 'make check -j'
            }
          }
        }
      }
    }
  }
}
