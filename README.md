# The featurizer application

## How to compile

  - update bioshell v.4 submodule:
  ```
  git submodule init
  git submodule update --rebase --remote
  ```
  
  - compile the app:
  ```
  cd featurizer
  cargo build --release
  ```
  
The featurizer app has been written in [rust](https://www.rust-lang.org/), you need to set up the toolchain
if you have never done that before. On Linux and macOS systems, this is done as follows:

``
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
``

To check if everything is fine, use the following in a terminal:
``
cargo -- version
``

Detailed description on installation is available [from rust documentation page](https://www.rust-lang.org/tools/install)
  