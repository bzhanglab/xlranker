import fs from "fs/promises";
import { exec } from "child_process";
import { promisify } from "util";

const execAsync = promisify(exec);

async function cleanAndBuild() {
  const dirPath = "../src/xlranker/static";

  try {
    await fs.rm(dirPath, { recursive: true, force: true });
    console.log(`Deleted ${dirPath}`);
  } catch (err) {
    console.error(`Error deleting ${dirPath}:`, err);
  }

  try {
    await fs.mkdir(dirPath, { recursive: true });
    console.log(`Created directory ${dirPath}`);

    await fs.writeFile(`${dirPath}/__init__.py`, "");
    console.log(`Created blank __init__.py`);
  } catch (err) {
    console.error("Error creating directory or file:", err);
    process.exit(1);
  }

  try {
    console.log("Running vite build...");
    const { stdout, stderr } = await execAsync("vite build");
    console.log(stdout);
    if (stderr) console.error(stderr);
  } catch (err) {
    console.error("Error running vite build:", err);
    process.exit(1);
  }
}

cleanAndBuild();
