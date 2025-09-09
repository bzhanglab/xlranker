<script lang="ts">
  import "./app.css";
  import Navbar from "$lib/components/Navbar.svelte";
  import { Label } from "$lib/components/ui/label";
  import { Input } from "$lib/components/ui/input/index.js";

  let peptide_file: File;

  async function check_network() {
    const name = peptide_file.name;
    document.getElementById("status").textContent = "Starting task...";
    document.getElementById("result").textContent = "";

    const response = await fetch("/start_task", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ name }),
    });

    const data = await response.json();
    currentTask = data.task_id;

    document.getElementById("status").textContent =
      "Task started. Checking progress...";
    pollTask();
  }
</script>

<Navbar />

<div class="container mx-auto p-4 prose prose-lg">
  <h1 class="text-center">XLRanker</h1>

  <h3>1. Peptide List</h3>
  <p>
    Upload your peptide sequences in a two column format separated by tabs or
    commas.
  </p>
  <div class="grid w-full max-w-sm items-center gap-1.5">
    <Label for="peptide-list">Peptide List</Label>
    <Input
      bind:value={peptide_file}
      id="peptide-list"
      type="file"
      oninput={check_network}
    />
  </div>
</div>
